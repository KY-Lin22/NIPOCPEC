classdef Hopper < DifferentialVariationalInequalities
    %% Hopper: 2D hopping robot with contact between foot and ground
    %   inspired by thowell's paper and source code: https://github.com/thowell/motion_planning/blob/main/models/hopper.jl
    %   1: lateral body position; 
    %   2; vertical body position;
    %   3: body orientation
    %   4: leg length
    properties
        mass (2, 1) double {mustBePositive, mustBeFinite} = [1; 0.1];
        inertia (2,1) double = [0.25; 0.025];
        mu = 1;% friction coefficient in ground
    end
    
    %% Constructor Method for Hopper   
    methods
        function plant = Hopper(timeStep, mass, inertia, mu)
            %Hopper
            %   Constructor of Class Hopper
            %           plant = Hopper(timeStep, mass, inertia, mu)
            %
            % Argument: 
            %           timeStep: (1 x 1) double          
            %           mass: (2 x 1) double
            %           inertia: (2 x 1) double 
            %           mu: (1 x 1) double
            %
            % Output: 
            %           plant: a object of Class Hopper 
            %% construct an object
            tau_Dim = 3;
            x_Dim = 8;
            p_Dim = 3;
            plant = plant@DifferentialVariationalInequalities(tau_Dim, x_Dim, p_Dim);
            
            % initialize dynamics variable limit
            plant.DynVarLimit.tau_Max = [50; 50; 100];
            plant.DynVarLimit.tau_Min = [-50; -50; -100]; 
            plant.DynVarLimit.x_Max = [2; 2; pi; 0.5; 10; 10; 5; 5]; 
            plant.DynVarLimit.x_Min = [0; 0; -pi; 0.1; -10; -10; -5; -5];   
            
            % initial other properties
            if nargin ~= 0
                if ~isempty(timeStep)
                    plant.timeStep = timeStep;
                end
                if ~isempty(mass)
                    plant.mass = mass;
                end
                if ~isempty(inertia)
                    plant.inertia = inertia;
                end 
                if ~isempty(mu)
                    plant.mu = mu;
                end
            end            
            %% set symbolic formulation about dynamics
            plant.codeOptimize = true;
            plant.setDynamics()                         
        end
        
    end
    %% Other Methods for Hopper
    methods
        % symbolic formulation about dynamics for codeGen
        setDynamics(plant)
        
        % Visualization Method
        plotConfiguration(plant, InitState, RefState)
        
        animateTrajectory(plant, simuTimeStep, InitState, tau, x, p)
        
        plotSimuResult(plant, simuTimeStep, InitState, tau, x, p)         
    end
end

