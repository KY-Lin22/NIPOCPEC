classdef CartPoleWithFriction < DifferentialVariationalInequalities
    %% CartPoleWithFriction: cart pole with friction bewteen cart and ground
    %   1: cart, 2: pole
    %
    properties
        mass (2, 1) double {mustBePositive, mustBeFinite} = [1; 0.1];
        linkLength (1, 1) double {mustBeNonnegative, mustBeFinite} = 1;
        cartLength (1, 1) double {mustBeNonnegative, mustBeFinite} = 0.4;
        
        cartHeight (1, 1) double {mustBeNonnegative, mustBeFinite} = 0.4;% just for visualization
    end
    %% Constructor Method for CartPoleWithFriction
    methods
        function plant = CartPoleWithFriction(timeStep, mass, linkLength, cartLength)
            %CartPoleWithFriction
            %   Constructor of Class CartPoleWithFriction
            %
            % Syntax: 
            %           plant = CartPoleWithFriction(timeStep, mass, linkLength, cartLength)
            %
            % Argument: 
            %           timeStep: (1 x 1) double    
            %           mass: (2 x 1) double
            %           linkLength: (1 x 1) double   
            %           cartLength: (1 x 1) double   
            %
            % Output: 
            %           plant: a object of Class CartPoleWithFriction            
            %
            %% construct an object
            tau_Dim = 1;
            x_Dim = 4;
            p_Dim = 1;
            plant = plant@DifferentialVariationalInequalities(tau_Dim, x_Dim, p_Dim);

            % initialize dynamics variable limit
            plant.DynVarLimit.tau_Max = 20; 
            plant.DynVarLimit.tau_Min = -20; 
            plant.DynVarLimit.x_Max = [5; 180/180*pi; 2; 1]; 
            plant.DynVarLimit.x_Min = [0; -180/180*pi; -2; -1];            
            
            % initial other properties
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
                if ~isempty(cartLength)
                    plant.cartLength = cartLength;
                end
            end
            %% set symbolic formulation about dynamics
            plant.codeOptimize = true;
            plant.setDynamics()            
        end
        
    end
    %% Other Methods for CartPoleWithFriction
    methods
        % symbolic formulation about dynamics for codeGen
        setDynamics(plant)
        
        % Visualization Method
        plotConfiguration(plant, InitState, RefState)
        
        animateTrajectory(plant, simuTimeStep, InitState, tau, x, p)
        
        plotSimuResult(plant, simuTimeStep, InitState, tau, x, p)        
        
    end
end

