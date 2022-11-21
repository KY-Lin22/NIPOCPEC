classdef FilippovDI_BVI < DifferentialVariationalInequalities
    %% FilippovDI_BVI: an example from Stewart's paper:
    %   Optimal control of systems with discontinuous differential equations
    %   
    
    properties
        
    end
    
    %% Constructor Method for FilippovDI_BVI
    methods
        function plant = FilippovDI_BVI(timeStep)
            %FilippovDI
            %   Constructor of Class FilippovDI_BVI
            %
            % Syntax: 
            %           plant = FilippovDI_BVI(timeStep)
            % Argument: 
            %           timeStep: (1 x 1) double 
            %
            % Output: 
            %           plant: a object of Class FilippovDI_BVI
            
            %% construct an object
            tau_Dim = 0;  
            x_Dim = 1;
            p_Dim = 1; % y
            plant = plant@DifferentialVariationalInequalities(tau_Dim, x_Dim, p_Dim);            
            
            % initialize dynamics variable limit
            plant.DynVarLimit.tau_Max = 10000*zeros(tau_Dim, 1);
            plant.DynVarLimit.tau_Min = -10000*zeros(tau_Dim, 1); 
            plant.DynVarLimit.x_Max = 10000; 
            plant.DynVarLimit.x_Min = -10000;
            
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
    
    %% Other Methods for FilippovDI_BVI
    methods
        % symbolic formulation about dynamics for codeGen
        setDynamics(plant)
        
        % Visualization Method
        plotSimuResult(plant, simuTimeStep, InitState, tau, x, p)
        
    end
end

