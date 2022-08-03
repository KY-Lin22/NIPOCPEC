classdef DifferentialVariationalInequalities < handle
    %% DifferentialVariationalInequalities: Superclass to construct dynamics system in the form of DVI
    %   We specify DVI as a dynamics version of BVI
    %   
    
    properties
        Dim; % struct, record dimension of variables
        
        tau; % (Dim.tau x 1) symbolic variable, control input
        x; % (Dim.x x 1) symbolic variable, system state
        p; % (Dim.p x 1) symbolic variable, equilibrium state
        
        DynVarLimit; % struct, limit for dynamics variables tau and x
    end
    
    properties
        codeOptimize; % logical, flag to determine whether optimizes code in matlabFunction
        computeInvM = false; % logical, flag for the computation of system dynamics dx = f(tau, x, p)
                             % (1)false: the formulation of state equation f(tau, x, p) is provided
                             %           hence it does not need to compute InvM explicitly 
                             % (2)true:  the formulation of state equation f(tau, x, p) is given by f = [dq;InvM*H]
                             %           hence it does need to compute InvM and H explicitly        
        f; % (Dim.x x 1) symbolic function, state equation
        M; % (0.5*Dim.x x 1) symbolic function, mass matrix
        H; % (0.5*Dim.x x 1) symbolic function, nonlinear function including centrifugal/Coriolis force, gravity and external force (control, contact force)
        K; % (Dim.p x 1) symbolic function, function in BVI(l, u, p, K(tau, x, p))
        l; % (Dim.p x 1) double, lower bound in BVI(l, u, p, K(tau, x, p))
        u; % (Dim.p x 1) double, upper bound in BVI(l, u, p, K(tau, x, p))
    end
    
    properties
        timeStep (1,1) double {mustBePositive, mustBeFinite} = 0.01;% discretization time step
        g(1,1) double = 9.8; % (1 x 1) double, gravity parameter 
    end
    
    %% Constructor Method for DifferentialVariationalInequalities     
    methods
        function plant = DifferentialVariationalInequalities(tau_Dim, x_Dim, p_Dim)
            %DifferentialVariationalInequalities
            %   Constructor of Class DifferentialVariationalInequalities
            %
            % Syntax: 
            %           plant = DifferentialVariationalInequalities(tau_Dim, x_Dim, p_Dim)
            % Argument:         
            %           tau_Dim: integer  
            %           x_Dim: integer
            %           p_Dim: interger
            %
            % Output: 
            %           plant: a object of Class DifferentialVariationalInequalities
            
            % check input argument
            mustBeInteger(tau_Dim)
            mustBeInteger(x_Dim)
            mustBeInteger(p_Dim)
            % assign properties value
            plant.Dim = struct('tau', tau_Dim,...
                               'x',   x_Dim,...
                               'p',   p_Dim);      
            % create symbolic variable to describe DVI
            plant.tau = sym('tau', [plant.Dim.tau, 1]);
            plant.x = sym('x', [plant.Dim.x, 1]);
            plant.p = sym('p', [plant.Dim.p, 1]);
            % initialize dynamics variable limit
            plant.DynVarLimit = struct('tau_Max', Inf(plant.Dim.tau, 1), 'tau_Min', -Inf(plant.Dim.tau, 1),...
                                      'x_Max',Inf(plant.Dim.x, 1), 'x_Min', -Inf(plant.Dim.x, 1));              
            % make dir and add to path
            mkdir('./autoGen_CodeFiles')
            addpath(genpath('./autoGen_CodeFiles'))
        end
        
    end
    
    %% Other Methods for DifferentialVariationalInequalities
    methods
        setDynVarLimit(plant, tau_Max, tau_Min, x_Max, x_Min) % set dynamics variable Limit
        
        % code generation and compututation of dynamics         
        codeGen(plant)
        
        f = computeStateEquation_Function(plant, tau, x, p)
        
        [ftau, fx, fp] = computeStateEquation_Jacobian(plant, tau, x, p)
        
        K = computeVI_Function(plant, tau, x, p)
    end
       
end

