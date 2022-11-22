classdef PeriodicDVI < DifferentialVariationalInequalities
    %% PeriodicDVI: combination of Holf oscillator and Filippov DI
    %   Detailed explanation goes here
    
    properties
        convergeSpeed1 = 50 % control the speed for the oscillator to converge to the limit cycle (rho_1)
        convergeSpeed2 = 50 % control the speed for the oscillator to converge to the limit cycle (rho_2)
        radiusLimitCycle = 1 % radius of the limit cycle
        periodLimitCycle = 0.6 % period of the limit cycle
        dutyFactor = 0.5 % duty factor determining the fraction of rho_2 > 0 in a whole period
        radiusSmoothVaries = 50 % ensure radiusLimitCycle varies smoothly across different half planes
        InitPhase = 0; % determine the gait pattern
        k_N = 0.5 % weight parameter in foot velocity function (normal and tangential)
        k_T = 0.25
    end
    %% Constructor Method for PeriodicDVI
    methods
        function plant = PeriodicDVI(timeStep, OscillatorParam, FootVelFuncParam)
            %PeriodicDVI
            %   Detailed explanation goes here
            %% construct an object
            num_Osci = 1;
            if~isempty(OscillatorParam.InitPhase)
                num_Osci = length(OscillatorParam.InitPhase);
            end        
            tau_Dim = 0 * num_Osci; % virtual control for foot trajectory function
            x_Dim = 4 * num_Osci; % x1 - x2: rho, oscillator state variable
                           % x3 - x4: foot position (normal and tangential)
            p_Dim = 2 * num_Osci; % p1 - p2: auxiliary variable for transforming Filippov DI into DVI (normal and tangential contact velocity)
                           
            plant = plant@DifferentialVariationalInequalities(tau_Dim, x_Dim, p_Dim); 
            
            % initialize dynamics variable limit
            plant.DynVarLimit.tau_Max = 10 * ones(tau_Dim, 1);
            plant.DynVarLimit.tau_Min = -10 * ones(tau_Dim, 1); 
            plant.DynVarLimit.x_Max = 100 * ones(x_Dim, 1); 
            plant.DynVarLimit.x_Min = -100 * ones(x_Dim, 1);       
            
            % initialize other properties
            if nargin ~= 0
                if ~isempty(timeStep)
                    plant.timeStep = timeStep;
                end
                % oscillator
                if ~isempty(OscillatorParam.convergeSpeed1)
                    plant.convergeSpeed1 = OscillatorParam.convergeSpeed1;
                end
                if ~isempty(OscillatorParam.convergeSpeed2)
                    plant.convergeSpeed2 = OscillatorParam.convergeSpeed2;
                end
                if ~isempty(OscillatorParam.radiusLimitCycle)
                    plant.radiusLimitCycle = OscillatorParam.radiusLimitCycle;
                end 
                if ~isempty(OscillatorParam.periodLimitCycle)
                    plant.periodLimitCycle = OscillatorParam.periodLimitCycle;
                end
                if ~isempty(OscillatorParam.dutyFactor)
                    plant.dutyFactor = OscillatorParam.dutyFactor;
                end
                if ~isempty(OscillatorParam.radiusSmoothVaries)
                    plant.radiusSmoothVaries = OscillatorParam.radiusSmoothVaries;
                end
                if ~isempty(OscillatorParam.InitPhase)
                    plant.InitPhase = OscillatorParam.InitPhase;
                end
                if~isempty(FootVelFuncParam.k_N)
                    plant.k_N = FootVelFuncParam.k_N;
                end
                if~isempty(FootVelFuncParam.k_T)
                    plant.k_T = FootVelFuncParam.k_T;
                end
            end
            %% set symbolic formulation about dynamics
            plant.codeOptimize = true;
            plant.setDynamics()              
        end
        
    end
    %% Other Methods for PeriodicDVI
    methods
        % symbolic formulation about dynamics for codeGen
        setDynamics(plant)
        
        % system simulation
        [state, phase, radius] = systemSimulation(plant, InitState, nStages, timeStep)
        
        animateSystemSimulationResult(plant, InitState, state, nStages, timeStep) 

    end    
    
end

