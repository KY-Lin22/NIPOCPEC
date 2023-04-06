classdef Hopper_PeriodicGait < DifferentialVariationalInequalities
    %% Hopper_PeriodicGait: 2D hopping robot with periodic gait and contact between foot and ground
    %
    properties
        frictionCoeff = 1 % friction coefficient in ground
        mass (2, 1) double {mustBePositive, mustBeFinite} = [1; 0.1];
        inertia (2, 1) double = [0.25; 0.025];       
    end
    
    properties
        convergeSpeed1 = 50 % control the speed for the oscillator to converge to the limit cycle (rho_1)
        convergeSpeed2 = 50 % control the speed for the oscillator to converge to the limit cycle (rho_2)
        radiusLimitCycle = 1 % radius of the limit cycle
        periodLimitCycle = 0.6 % period of the limit cycle
        dutyFactor = 0.5 % duty factor determining the fraction of rho_2 > 0 in a whole period
        radiusSmoothVaries = 50 % ensure radiusLimitCycle varies smoothly across different half planes
    end
    
    %% Constructor Method for Hopper_PeriodicGait     
    methods
        function plant = Hopper_PeriodicGait(timeStep, frictionCoeff, BodyParam, OscillatorParam)
            %Hopper_PeriodicGait
            %   Detailed explanation goes here
            %% construct an object
            tau_Dim = 3;  % tau1: body moment taub
                          % tau2: leg force taul
                          % tau3: tangential friction
            x_Dim = 8 + 2 + 2; % x1: lateral body position; 
                               % x2; vertical body position;
                               % x3: body orientation
                               % x4: leg length
                               % x5 - x8:  velocity of x1 - x4
                               % x9 - x10: rho, oscillator state variable
                               % x11 - x12: foot position (normal and tangential)
            p_Dim = 3 + 2; % p1: normal contact force
                           % p2 - p3: auxiliary variable for tangential contact velocity p2 - p3 = vel_T;
                           % p4 - p5: auxiliary variable for transforming Filippov DI into DVI (normal and tangential contact velocity)
                           
            plant = plant@DifferentialVariationalInequalities(tau_Dim, x_Dim, p_Dim);
            
            % initialize dynamics variable limit
            plant.DynVarLimit.tau_Max = [50; 50; 100];
            plant.DynVarLimit.tau_Min = [-50; -50; -100]; 
            plant.DynVarLimit.x_Max = [2; 2; pi; 0.5; 10; 10; 5; 5;...
                100; 100; 100; 100]; 
            plant.DynVarLimit.x_Min = [0; 0; -pi; 0.1; -10; -10; -5; -5;...
                -100; -100; -100; -100]; 
            
            % initialize other properties
            if nargin ~= 0
                if ~isempty(timeStep)
                    plant.timeStep = timeStep;
                end
                if ~isempty(frictionCoeff)
                    plant.frictionCoeff = frictionCoeff;
                end
                % body
                if ~isempty(BodyParam.mass)
                    plant.mass = BodyParam.mass;
                end
                if ~isempty(BodyParam.inertia)
                    plant.inertia = BodyParam.inertia;
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
            end
            
            %% set symbolic formulation about dynamics
            plant.codeOptimize = true;
            plant.setDynamics()             
            
        end
        
    end
    %% Other Methods for Hopper_PeriodicGait
    methods
        % symbolic formulation about dynamics for codeGen
        setDynamics(plant)
        
        % Visualization Method
        plotConfiguration(plant, InitState, RefState)
        
        animateTrajectory(plant, simuTimeStep, InitState, tau, x, p)
        
        plotSimuResult(plant, simuTimeStep, InitState, tau, x, p)          
    end
    
end

