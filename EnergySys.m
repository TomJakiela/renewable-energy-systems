classdef EnergySys
    %ENERGYSYS Encapsulates power system operation
    %   Properties are expressed in kWh per time step
    %
    %   (c) Thomas Jakiela 2020 - 2021
    %   Available at https://github.com/TomJakiela/renewable-energy-systems
    %
    
    properties
        dailyCap        % in kWh total capacity e.g. battery size
        p2gCap          % in kWh per time step (power * sampling period/hr)
        p2gEta          % efficiency to gas energy content, not round trip
        g2pEta          % efficiency from gas to electricity
        
        dailykWh = 0    % Initial short-term storage e.g. battery charge
        seasonalkWh = 0 % Initial long-term storage e.g. CNG
        wastekWh = 0    % Curtailment
        
        deltaD = 0      % change in daily storage       % used during
        deltaS = 0      % change in seasonal storage    % simulation
        deltaC = 0      % change in curtailment         %

    end
    
    properties(Constant)
        % These properties are not inherent to the class but depend on the
        % locality and time resolution of the dataset. *Change as needed
        summer = 3:10;  % months for prioritizing PtG       *
        perHr = 4;      % Time samples per hour             *
    end

    methods
        function obj = EnergySys( ...
                dailykWhCap, power2gasCapacity, power2gasEta, gas2powerEta)
            %INIT Construct an instance of this class
            %
            %   dailykWhCap             Short-term storage (battery) size
            %   power2gasCapacity       PtG size in e-kWh per time step
            %   power2gasEta            PtG efficiency in kWh per e-kWh
            %   gas2powerEta            GtP efficiency in e-kWh per kWh
            
            %   Returns a struct holding the parameters
            %
            obj.dailyCap = dailykWhCap;
            obj.p2gCap = power2gasCapacity / obj.perHr;
            obj.p2gEta = power2gasEta;
            obj.g2pEta = gas2powerEta;
        end
        
        function attrs = getSums(obj)
            attrs = [obj.dailykWh obj.seasonalkWh obj.wastekWh];
        end
        
        function attrs = getDeltas(obj)
            attrs = [obj.deltaD obj.deltaS obj.deltaC];
        end
        
        function [obj, remainder] = power2gas(obj, powerInput, ctrl)
            %POWER2GAS Increments long term storage according to args
            %
            %   powerInput      net power available at the time step
            %   ctrl            signal to run p2g from battery or not
            %   
            %   remainder       excess or deficit after allocating energy 
            %                   this function does not directly withdraw
            %                   from battery, remainder is handled
            %                   externally
            %
            %   ctrl will be passed in if external signal determines there
            %   is a net surplus in the forecast. Function will not add
            %   load without surplus power input or available battery, the
            %   deficit will pass through to gas2power. ctrl is a range of
            %   values from 0 to 1
            %
            %   A more accurate model will use longer term forecasting to
            %   limit intermittent start and stop cycles, however this
            %   appears to be infrequent in practical AY simulations.
            
            minB = (1/10); % Minimum battery reserve to run power2gas
            if ctrl > 0
                if obj.dailykWh > (obj.dailyCap * minB)
                    p2g = min(obj.dailykWh - (obj.dailyCap * minB), ...
                        obj.p2gCap);
                else
                   p2g = 0;
                end
            else
                p2g = 0;
                % Simply no gas storage if no external command
            end
            remainder = powerInput - p2g;
            obj.deltaS = p2g;
            obj.seasonalkWh = obj.seasonalkWh + p2g * obj.p2gEta;
        end
        
        function [obj, remainder] = battery(obj, netEnergyDiff)
            %BATTERY Charges and discharges the short term storage.
            %
            %   Assumes ideal efficiency. Charges if positive input,
            %   discharges if negative input. 
            %
            if netEnergyDiff > 0
                deltaCharge = min(netEnergyDiff, obj.dailyCap - obj.dailykWh);
            else
                deltaCharge = max(-1 * obj.dailykWh, netEnergyDiff);
            end
            remainder = netEnergyDiff - deltaCharge;
            obj.deltaD = deltaCharge;
            obj.dailykWh = obj.dailykWh + deltaCharge;
        end
        
        function [fcast] = forecast(obj, inputs)
            %FORECAST Enables forecasting for control of power2gas
            %
            %   This operates on the length of simulation data to determine
            %   the net energy flow over the next 24 hours. Entire year
            %   forecast is not known in practice but is precomputed to
            %   simplify the simulation.  Noise or bias can be added to
            %   simulate meteorological uncertainty.
            %
            %   inputs = [gen, load] two column vectors of kWh deltas
            %   fcast = Matrix, row vectors are the length of the forecast
            %
            if obj.p2gCap > 0
                stepsAhead = floor(obj.perHr * 24 * 1);
                gen = [inputs(:, 1); zeros(stepsAhead, 1)];
                load = [inputs(:, 2); zeros(stepsAhead, 1)];
                fcast = zeros(size(inputs, 1), stepsAhead);
                for dT = 1:size(inputs, 1)
                    timesteps = dT:(dT + stepsAhead - 1);
                    fcast(dT, :) = cumsum(( ...
                        gen(timesteps) - load(timesteps)))';
                end
            else
                % no seasonal storage (PtG), skip forecasting
                fcast = inputs * 0;
            end   
        end
        
        function obj = step(obj, inputs, fcast)
            %STEP Iterates a single time step through the data
            %
            %   inputs      [gen load] in kWh per time step
            %   fcast       net kWh forecast looking forward 24 hours
            %
            %   returns none, modifies storage and delta attributes
            %
            %   power2gas either adds load or passes deficit through.
            %   battery attempts to cover remaining load or absorb surplus
            %   remaining surplus is discarded as curtailment
            %   remaining e-load is withdrawn from gas at rate g2pEta
            %   heat load always comes from gas storage
            %
            gen = inputs(1); eload = inputs(2);
            heatload = inputs(3); month = inputs(4);
            
            mismatch = gen - eload;
            if any(month == obj.summer)
                fullp2g = (... % This generates a cumulitive sum of load
                    1:length(fcast(:))) * obj.p2gCap;
                seasonCtrl = max( ...
                    all((fcast + obj.dailykWh - fullp2g ) > 0), ...
                    all((fcast + obj.dailykWh - fullp2g*.25) > 0) * .25);
            else
                seasonCtrl = 0;
            end
            [obj, rem] = obj.power2gas(mismatch, seasonCtrl);
            [obj, rem] = obj.battery(rem);
            
            if rem < 0
                obj.deltaS = rem;
                obj.seasonalkWh = obj.seasonalkWh + ...
                    rem * obj.g2pEta; % gas to electric efficiency
                % Slow startup and operation is neglected. Battery storage
                % and forecasting would be expected to provide sufficient
                % notice for startup in practice
            else
                obj.wastekWh = obj.wastekWh + rem;
                obj.deltaC = rem;
            end
             % Remaining additions to this incarnation
            obj.seasonalkWh = obj.seasonalkWh - heatload;
        end
        
        function [kWhSum, kWhDelta] = run(obj, inputData)
            %RUN Executes simulation of power system
            %
            %   inputData = [gen, e-load, month indicator, heat-load]
            %
            %   kWhSum   =  [battery storage, seasonal storage]
            %   kWhDelta ~~  diff([battery, seasonal, curtailment])
            %
            kWhSum = zeros(length(inputData), length(obj.getSums()));
            kWhDelta = zeros(length(inputData), length(obj.getDeltas()));
            fcast = obj.forecast(inputData);
            for x = 1:size(inputData, 1)
                % x indexes the input rows
                obj = obj.step(inputData(x, :)', fcast(x, :));
                kWhSum(x, :) = obj.getSums();
                kWhDelta(x, :) = obj.getDeltas();
            end
        end
        
    end % end methods
end % end class
