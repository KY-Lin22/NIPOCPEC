classdef FilippovDI_control < DifferentialVariationalInequalities
    %% FilippovDI_control: an example extended from Stewart's paper:
    %   Optimal control of systems with discontinuous differential equations
    %   
    
    properties
        
    end
    
    %% Constructor Method for FilippovDI_control
    methods
        function plant = FilippovDI_control(timeStep)
            %FilippovDI_control
            %   Constructor of Class FilippovDI_control
            %
            % Syntax: 
            %           plant = FilippovDI_control(timeStep)
            % Argument: 
            %           timeStep: (1 x 1) double 
            %
            % Output: 
            %           plant: a object of Class FilippovDI_control
            
            %% construct an object
            tau_Dim = 2; % tau1: y
                         % tau2: control
            x_Dim = 1;
            p_Dim = 2; % p1: lambda0
                       % p2: lambda1
            plant = plant@DifferentialVariationalInequalities(tau_Dim, x_Dim, p_Dim);            
            
            % initialize dynamics variable limit
            plant.DynVarLimit.tau_Max = [10000; 5];
            plant.DynVarLimit.tau_Min = [-10000; 0]; 
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
    
    %% Other Methods for FilippovDI_control
    methods
        % symbolic formulation about dynamics for codeGen
        setDynamics(plant)
        
        % Visualization Method
        plotSimuResult(plant, simuTimeStep, InitState, tau, x, p)
        
    end
end