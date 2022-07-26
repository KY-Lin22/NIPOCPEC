classdef Acrobot < DifferentialVariationalInequalities
    %% Acrobot: acrobot with contact due to the joint limit in joint 2
    %   Detailed explanation goes here
    
    properties
        mass (2,1) double {mustBePositive, mustBeFinite} = [1; 1]; 
        linkLength (2,1) double {mustBeNonnegative, mustBeFinite} = [1; 1]; 
        linkCenter (2,1) double = [0.5; 0.5]; 
        inertia (2,1) double = [1; 1]; 
        jointFriction (2, 1) double = [0.1; 0.1];
        q2_min (1, 1) double {mustBeFinite} = -pi;
        q2_max (1, 1) double {mustBeFinite} = pi; 
    end
    %% Constructor Method for Acrobot  
    methods
        function plant = Acrobot(timeStep, mass, linkLength, linkCenter, inertia, q2_min, q2_max)
            %Acrobot
            %   Constructor of Class Acrobot
            %
            % Syntax: 
            %           plant = Acrobot(timeStep, mass, linkLength, linkCenter, inertia, q2_min, q2_max)
            %
            % Argument: 
            %           timeStep: (1 x 1) double          
            %           mass: (2 x 1) double
            %           linkLength: (2 x 1) double  
            %           linkCenter: (2 x 1) double  
            %           inertia: (2 x 1) double 
            %           q2_min: (1 x 1) double
            %           q2_max: (1 x 1) double
            %
            % Output: 
            %           plant: a object of Class Acrobot
            %% construct an object
            tau_Dim = 1;
            x_Dim = 4;
            p_Dim = 2;                        
            plant = plant@DifferentialVariationalInequalities(tau_Dim, x_Dim, p_Dim);
            
            % initialize dynamics variable limit
            plant.DynVarLimit.tau_Max = 50;
            plant.DynVarLimit.tau_Min = -50; 
            plant.DynVarLimit.x_Max = [2*pi; 1/2*pi; 5; 5]; 
            plant.DynVarLimit.x_Min = [-2*pi; -1/2*pi; -5; -5];            
            
            % initialize other properties
            if nargin ~= 0
                if ~isempty(timeStep)
                    plant.timeStep = timeStep;
                end
                if ~isempty(mass)
                    plant.mass = mass;
                end
                if ~isempty(linkLength)
                    plant.linkLength = linkLength;
                end
                if ~isempty(linkCenter)
                    plant.linkCenter = linkCenter;
                end
                if ~isempty(inertia)
                    plant.inertia = inertia;
                end
                if ~isempty(q2_min)
                    plant.q2_min = q2_min;
                    plant.DynVarLimit.x_Min(2) = q2_min;
                end
                if ~isempty(q2_max)
                    plant.q2_max = q2_max;
                    plant.DynVarLimit.x_Max(2) = q2_max;
                end
            end
            %% set symbolic formulation about dynamics
            plant.codeOptimize = true;
            plant.setDynamics()              
            
        end
       
    end
    %% Other Methods for Acrobot
    methods
        % symbolic formulation about dynamics for codeGen
        setDynamics(plant)
        
        % Visualization Method
        plotConfiguration(plant, InitState, RefState)
        
        animateTrajectory(plant, simuTimeStep, InitState, tau, x, p)
        
        plotSimuResult(plant, simuTimeStep, InitState, tau, x, p)
    end
end

