classdef ThreeCarts < DifferentialVariationalInequalities
    %% ThreeCarts: a three carts system with contact bewteen carts
    %   Detailed explanation goes here
    
    properties
        mass (3,1) double {mustBePositive, mustBeFinite} = [1; 1; 1]; 
        cartLength (3,1) double {mustBeNonnegative, mustBeFinite} = [2; 2; 2]
        cartHeight (3,1) double {mustBeNonnegative, mustBeFinite} = [2; 2; 2] % just for visualization
        viscousDamping (1, 1) {mustBePositive, mustBeFinite} = 2;
    end
    %% Constructor Method for ThreeCarts
    methods
        function plant = ThreeCarts(timeStep, mass, cartLength, viscousDamping)
            %ThreeCarts
            %   Constructor of Class ThreeCarts
            %
            % Syntax: 
            %           plant = ThreeCarts(timeStep, mass, cartLength, viscousDamping)
            %
            % Argument: 
            %           timeStep: (1 x 1) double    
            %           mass: (3 x 1) double
            %           cartLength: (3 x 1) double              
            %           viscousDamping: (1 x 1) double 
            %
            % Output: 
            %           plant: a object of Class ThreeCarts 
            %
            %   Detailed explanation goes here
            %% construct an object
            tau_Dim = 3;
            x_Dim = 6;
            p_Dim = 2;
            plant = plant@DifferentialVariationalInequalities(tau_Dim, x_Dim, p_Dim);            
            
            % initialize dynamics variable limit
            plant.DynVarLimit.tau_Max = [20; 20; 20]; 
            plant.DynVarLimit.tau_Min = [-20; -20; -20]; 
            plant.DynVarLimit.x_Max = [10; 15; 15; 5; 5; 5]; 
            plant.DynVarLimit.x_Min = [-15; -15; -10; -5; -5; -5];            
            
            % initialize other properties
            if nargin ~= 0
                if ~isempty(timeStep)
                    plant.timeStep = timeStep;
                end
                if ~isempty(mass)
                    plant.mass = mass;
                end
                if ~isempty(cartLength)
                    plant.cartLength = cartLength;    
                end
                if ~isempty(viscousDamping)
                    plant.viscousDamping = viscousDamping;
                end            
            end             
            %% set symbolic formulation about dynamics
            plant.codeOptimize = true;
            plant.setDynamics()             
            
        end
    end
    
    %% Other Methods for ThreeCarts
    methods
        % symbolic formulation about dynamics for codeGen
        setDynamics(plant)
        
        % Visualization Method
        plotConfiguration(plant, InitState, RefState)
        
        animateTrajectory(plant, simuTimeStep, InitState, tau, x, p)
        
        plotSimuResult(plant, simuTimeStep, InitState, tau, x, p)         
    end    
    
end

