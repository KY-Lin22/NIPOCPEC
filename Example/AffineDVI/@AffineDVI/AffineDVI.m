classdef AffineDVI < DifferentialVariationalInequalities
    %% AffineDVI: an example generalized from Vieira's 19TAC paper:
    %  Quadratic Optimal Control of Linear Complementarity Systems: First order necessary conditions and numerical analysis
    
    properties
        
    end
    %% Constructor Method for AffineDVI
    methods
        function plant = AffineDVI(timeStep)
            %AffineDVI
            %   Constructor of Class AffineDVI
            %
            % Syntax: 
            %           plant = AffineDVI(timeStep)
            % Argument: 
            %           timeStep: (1 x 1) double 
            %
            % Output: 
            %           plant: a object of Class AffineDVI
            %
            %% construct an object
            tau_Dim = 1;
            x_Dim = 2;
            p_Dim = 1;
            plant = plant@DifferentialVariationalInequalities(tau_Dim, x_Dim, p_Dim);
            
            % initialize dynamics variable limit
            plant.DynVarLimit.tau_Max = 20;
            plant.DynVarLimit.tau_Min = -20; 
            plant.DynVarLimit.x_Max = [5; 5]; 
            plant.DynVarLimit.x_Min = [-5; -5];             
            
           % initialize timeStep
            if nargin ~= 0
                if ~isempty(timeStep)
                    plant.timeStep = timeStep;
                end          
            end            
            
            %% set symbolic formulation about dynamics
            plant.codeOptimize = true;
            plant.setDynamics()
            
        end

    end
    
    %% Other Methods for AffineDVI
    methods
        % symbolic formulation about dynamics for codeGen
        setDynamics(plant)
        
        % Visualization Method
        plotConfiguration(plant, InitState, RefState)
        
        animateTrajectory(plant, simuTimeStep, InitState, tau, x, p)
        
        plotSimuResult(plant, simuTimeStep, InitState, tau, x, p)
        
    end
    
end

