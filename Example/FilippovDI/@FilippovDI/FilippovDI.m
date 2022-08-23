classdef FilippovDI < DifferentialVariationalInequalities
    %% FilippovDI: an example from Stewart's paper:
    %   Optimal control of systems with discontinuous differential equations
    %   
    
    properties
        
    end
    
    %% Constructor Method for FilippovDI
    methods
        function plant = FilippovDI(timeStep)
            %FilippovDI
            %   Constructor of Class FilippovDI
            %
            % Syntax: 
            %           plant = FilippovDI(timeStep)
            % Argument: 
            %           timeStep: (1 x 1) double 
            %
            % Output: 
            %           plant: a object of Class FilippovDI
            
            %% construct an object
            tau_Dim = 1; % tau1: y
            x_Dim = 1;
            p_Dim = 2; % p1: lambda0
                       % p2: lambda1
            plant = plant@DifferentialVariationalInequalities(tau_Dim, x_Dim, p_Dim);            
            
            % initialize dynamics variable limit
            plant.DynVarLimit.tau_Max = 10000;
            plant.DynVarLimit.tau_Min = -10000; 
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
    
    %% Other Methods for FilippovDI
    methods
        % symbolic formulation about dynamics for codeGen
        setDynamics(plant)
        
        % Visualization Method
        plotSimuResult(plant, simuTimeStep, InitState, tau, x, p)
        
    end
end

