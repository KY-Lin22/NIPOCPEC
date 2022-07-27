classdef OCPEC_Perturbed < handle
    %% OCPEC_Perturbed: Superclass to construct perturbed OCPEC
    %   We provide three kinds of method to reformulate the original OCPEC as a perturbed OCPEC
    %   1.SmoothingEquation
    %     We reformulate VI as the smoothing approximation of nonsmooth normal equation
    %     here we choose CHKS smoothing function with perturbed parameter s
    %   2.Reg_NCPs
    %     We reformulate VI as the regularized nonlinear complementary problem 
    %   3.Reg_Scholtes
    %     We reformulate VI as the Scholtes regularization reformulation, see paper:
    %    'Convergence Properties of a Regularization Scheme for Mathematical Programs with Complementarity Constraints'
    %    2001, SIAM Journal on Optimization    
    %
    
    properties
        plant; % object, a object of DifferentialVariationalInequalities 
        timeStep (1,1) double {mustBePositive, mustBeFinite} = 0.01; % discretization time step 
        nStages {mustBeInteger}; % number of discretized stages
        InitState; % double, given initial state x0
        StageCost; % struct, with field 'xRef', 'tauRef', 'xWeight', 'tauWeight'
        TerminalCost; % struct, with field 'xRef', 'tauRef', 'xWeight', 'tauWeight'
        VI_mode char {mustBeMember(VI_mode, {'SmoothingEquation', 'Reg_NCPs', 'Reg_Scholtes'})} = 'Reg_Scholtes'; % types of VI perturbed reformulation
    end
    
    properties
        Dim; % struct, record dimension of variables
        
        tau; % (Dim.tau x 1) sym var, control input
        x; % (Dim.x x 1) sym var, system state
        p; % (Dim.p x 1) sym var, equilibrium state        
        w; % (Dim.w x 1) sym var, auxiliary variable for function K in BVI
        
        sigma; % (Dim.sigma x 1) sym var, dual variable for inequality constraints G        
        eta; % (Dim.eta x 1) sym var, dual variable for equality constraints C
        lambda; % (Dim.lambda x 1) sym var, dual variable for discretized system dynamics F
        gamma; % (Dim.gamma x 1) sym var, dual variable for perturbed equilibrium dynamics PHI 
    end  
    
    properties
        codeOptimize; % logical, flag to determine whether optimizes code in matlabFunction                 
            
        symVar; % struct, record auxiliary sym var used in sym fun's formulation and code generation 
        LS; % sym fun, stage cost function LS(tau, x, p, auxivar)
        LT; % sym fun, terminal cost function LT(tau, x, p, auxivar)
        G; % sym fun, inequality constraints G (tau, x, p) >= 0          
        C; % sym fun, equality constraints C(tau, x, p, auxivar) = 0         
        PHI; % sym fun, perturbed equilibrium dynamics PHI(p, auxivar, s) with perturbed parameter s >= 0, 
        
        LFRP; % sym fun, stage cost function LFRP(tau, x, p, auxivar) used in Feasibility Restoration Phase
        FRP; % struct, with field 'ZRef' and 'ZWeight', used to record ref and weight matrix for L_FRP           
    end

    
%% Constructor Method for OCPEC_Perturbed    
    methods
        function OCPEC = OCPEC_Perturbed(plant, timeStep, nStages, InitState, StageCost, TerminalCost, VI_mode)
            %OCPEC_Perturbed
            %   Constructor of Class OCPEC_Perturbed
            %
            % Syntax: 
            %           OCPEC = OCPEC_Perturbed(plant, timeStep, nStages, InitState, StageCost, TerminalCost, VI_mode)
            % Argument: 
            %           plant: obj
            %           timeStep: (1 x 1) double
            %           nStages: integer
            %           InitState: (Dim.x x 1) double, 
            %           StageCost: struct, with field 'xRef'(Dim.x x nStages), 'tauRef'(Dim.tau x nStages), 'xWeight'(Dim.x x 1), 'tauWeight'(Dim.tau x 1)
            %           TerminalCost: struct, with field 'xRef'(Dim.x x 1), 'tauRef'(Dim.tau x 1), 'xWeight'(Dim.x x 1), 'tauWeight'(Dim.tau x 1)
            %           VI_mode: char, can be assigned as 'SmoothingEquation', 'Reg_NCPs', 'Reg_Scholtes'
            %
            % Output: 
            %           OCPEC: a object of Class OCPEC_Perturbed
            %
            %% check input      
            if ~all(size(InitState) == [plant.Dim.x, 1])
                error('please specify the correct InitState')
            end
            % check stageCost
            if ~all(size(StageCost.xRef) == [plant.Dim.x, nStages])
                error('please specify the correct StageCost.xRef')
            end
            if ~all(size(StageCost.tauRef) == [plant.Dim.tau, nStages])
                error('please specify the correct StageCost.tauRef')
            end   
            if ~all(size(StageCost.xWeight) == [plant.Dim.x, 1])
                error('please specify the correct StageCost.xWeight')
            end
            if ~all(size(StageCost.tauWeight) == [plant.Dim.tau, 1])
                error('please specify the correct StageCost.tauWeight')
            end
            % check terminalCost
            if ~all(size(TerminalCost.xRef) == [plant.Dim.x, 1])
                error('please specify the correct TerminalCost.xRef')
            end
            if ~all(size(TerminalCost.tauRef) == [plant.Dim.tau, 1])
                error('please specify the correct TerminalCost.tauRef')
            end   
            if ~all(size(TerminalCost.xWeight) == [plant.Dim.x, 1])
                error('please specify the correct TerminalCost.xWeight')
            end
            if ~all(size(TerminalCost.tauWeight) == [plant.Dim.tau, 1])
                error('please specify the correct TerminalCost.tauWeight')
            end
            
            %% initialize properties: basic         
            OCPEC.plant = plant;
            if ~isempty(timeStep)
                OCPEC.timeStep = timeStep;
            else
                OCPEC.timeStep = plant.timeStep;
            end         
            OCPEC.nStages = nStages;
            OCPEC.InitState = InitState;           
            OCPEC.StageCost = StageCost;               
            OCPEC.TerminalCost = TerminalCost;   
            OCPEC.VI_mode = VI_mode;            
            
            %% initialize properties: primal variable
            % initialize basic variable tau, x, p
            OCPEC.Dim = struct('tau', plant.Dim.tau,...
                              'x', plant.Dim.x,...
                              'p', plant.Dim.p);
            OCPEC.tau = plant.tau;
            OCPEC.x   = plant.x;
            OCPEC.p   = plant.p;             
            % initialize auxiliary variable w
            NCP_num = 0;
            for i = 1 : OCPEC.Dim.p
                if (OCPEC.plant.l(i) == 0) && (OCPEC.plant.u(i) == Inf)
                    NCP_num = NCP_num + 1;
                end
            end
            BVI_num = OCPEC.Dim.p - NCP_num;           
            if (strcmp(OCPEC.VI_mode, 'SmoothingEquation')) || (strcmp(OCPEC.VI_mode, 'Reg_Scholtes'))
                OCPEC.Dim.w = OCPEC.Dim.p;
            elseif strcmp(OCPEC.VI_mode, 'Reg_NCPs')
                OCPEC.Dim.w = NCP_num + 2 * BVI_num;
            end
            OCPEC.w = sym('w', [OCPEC.Dim.w, 1]);
            OCPEC.Dim.Z = OCPEC.Dim.tau + OCPEC.Dim.x + OCPEC.Dim.p + OCPEC.Dim.w;
            
            %% initialize properties: cost function
            OCPEC.codeOptimize = plant.codeOptimize;             
            OCPEC.symVar = struct();            
            
            % sym var for ref and weight matrix used to construct cost function L
            OCPEC.symVar.xRef = sym('xRef', [OCPEC.Dim.x, 1]);
            OCPEC.symVar.tauRef = sym('tauRef', [OCPEC.Dim.tau, 1]);
            OCPEC.symVar.xWeight = sym('xWeight', [OCPEC.Dim.x, 1]); % diagonal elements
            OCPEC.symVar.tauWeight = sym('tauWeight', [OCPEC.Dim.tau, 1]); % diagonal elements
            
            % stage cost function LS
            pWeight = 0.001 * eye(OCPEC.Dim.p);
            wWeight = 0.001 * eye(OCPEC.Dim.w);
            LS_x = 0.5 * (OCPEC.x - OCPEC.symVar.xRef)' * diag(OCPEC.symVar.xWeight) * (OCPEC.x - OCPEC.symVar.xRef);
            LS_tau = 0.5 * (OCPEC.tau - OCPEC.symVar.tauRef)' * diag(OCPEC.symVar.tauWeight) * (OCPEC.tau - OCPEC.symVar.tauRef);
            LS_p = 0.5 * OCPEC.p' * pWeight * OCPEC.p;
            LS_w = 0.5 * OCPEC.w' * wWeight * OCPEC.w;
            OCPEC.LS = (LS_x + LS_tau + LS_p + LS_w) * OCPEC.timeStep;
            
            % terminal cost function LT
            LT_x = 0.5 * (OCPEC.x - OCPEC.symVar.xRef)' * diag(OCPEC.symVar.xWeight) * (OCPEC.x - OCPEC.symVar.xRef);
            LT_tau = 0.5 * (OCPEC.tau - OCPEC.symVar.tauRef)' * diag(OCPEC.symVar.tauWeight) * (OCPEC.tau - OCPEC.symVar.tauRef);
            OCPEC.LT = LT_x + LT_tau;
            
            %% initialize properties: constraints function G and C
            % inequality constraints G
            DynVarLimit = plant.DynVarLimit;
            OCPEC.G = [DynVarLimit.tau_Max - OCPEC.tau;...
                OCPEC.tau - DynVarLimit.tau_Min;...
                DynVarLimit.x_Max - OCPEC.x;...
                OCPEC.x - DynVarLimit.x_Min];
                   
            % equality constraint C
            if (strcmp(OCPEC.VI_mode, 'SmoothingEquation')) || (strcmp(OCPEC.VI_mode, 'Reg_Scholtes'))
                % define auxiliary variable w = K
                auxiVarEquation = OCPEC.w - OCPEC.plant.K;
            elseif strcmp(OCPEC.VI_mode, 'Reg_NCPs')
                % define auxiliary variable NCP: w1 = K(i), BVI: w1-w2 = K(i)
                auxiVarEquation = sym('auxiVarEquation', [OCPEC.Dim.p, 1]);
                w_Counter = 0;
                for i = 1 : OCPEC.Dim.p
                    if (OCPEC.plant.l(i) == 0) && (OCPEC.plant.u(i) == Inf)
                        % nonlinear complementary problem
                        auxiVarEquation(i, 1) = OCPEC.w(w_Counter + 1, 1) - OCPEC.plant.K(i);
                        w_Counter = w_Counter + 1;
                    else
                        % box constraint variation inequality
                        auxiVarEquation(i, 1) = OCPEC.w(w_Counter + 1, 1) - OCPEC.w(w_Counter + 2, 1) - OCPEC.plant.K(i);
                        w_Counter = w_Counter + 2;
                    end
                end
            end
            OCPEC.C = auxiVarEquation;
            
            %% initialize properties: constraints function PHI
            OCPEC.symVar.s = sym('s', [1, 1]); % perturbed parameter used in PHI
            switch OCPEC.VI_mode
                case 'SmoothingEquation'
                    OCPEC.PHI = sym('PHI', [OCPEC.Dim.p, 1]);
                    for i = 1 : OCPEC.Dim.p
                        if (OCPEC.plant.l(i) == 0) && (OCPEC.plant.u(i) == Inf)
                            % nonlinear complementary problem
                            pw = OCPEC.p(i) - OCPEC.w(i);
                            h_i = 1/2 * (sqrt(pw^2 + 4 * OCPEC.symVar.s^2) + pw);
                        else
                            % box constraint variation inequality
                            lpw = OCPEC.plant.l(i) - (OCPEC.p(i) - OCPEC.w(i));
                            upw = OCPEC.plant.u(i) - (OCPEC.p(i) - OCPEC.w(i));
                            h_i = 1/2*(OCPEC.plant.l(i) + sqrt(lpw^2 + 4*OCPEC.symVar.s^2)) ...
                                + 1/2*(OCPEC.plant.u(i) - sqrt(upw^2 + 4*OCPEC.symVar.s^2));
                        end
                        OCPEC.PHI(i, 1) = OCPEC.p(i) - h_i;
                    end
                case 'Reg_NCPs'
                    OCPEC.PHI = sym('PHI', [3 * NCP_num + 6 * BVI_num, 1]);
                    w_Counter = 0;
                    PHI_Counter = 0;
                    for i = 1 : OCPEC.Dim.p
                        if (OCPEC.plant.l(i) == 0) && (OCPEC.plant.u(i) == Inf)
                            % nonlinear complementary problem
                            OCPEC.PHI(PHI_Counter + 1 : PHI_Counter + 3, 1) =...
                                [OCPEC.p(i);...
                                OCPEC.w(w_Counter + 1, 1);...
                                OCPEC.symVar.s - OCPEC.p(i) * OCPEC.w(w_Counter + 1, 1)];
                            w_Counter = w_Counter + 1;
                            PHI_Counter = PHI_Counter + 3;
                        else
                            OCPEC.PHI(PHI_Counter + 1 : PHI_Counter + 6, 1) =...
                                [OCPEC.p(i) - OCPEC.plant.l(i);...
                                OCPEC.w(w_Counter + 1, 1);...
                                OCPEC.symVar.s - (OCPEC.p(i) - OCPEC.plant.l(i)) * OCPEC.w(w_Counter + 1, 1);...
                                OCPEC.plant.u(i) - OCPEC.p(i);...
                                OCPEC.w(w_Counter + 2, 1);...
                                OCPEC.symVar.s - (OCPEC.plant.u(i) - OCPEC.p(i)) * OCPEC.w(w_Counter + 2, 1)];
                            w_Counter = w_Counter + 2;
                            PHI_Counter = PHI_Counter + 6;
                        end
                    end                   
                case 'Reg_Scholtes'
                    OCPEC.PHI = sym('PHI', [3 * NCP_num + 4 * BVI_num, 1]);
                    PHI_Counter = 0;
                    for i = 1 : OCPEC.Dim.p
                        if (OCPEC.plant.l(i) == 0) && (OCPEC.plant.u(i) == Inf)
                            % nonlinear complementary problem
                            OCPEC.PHI(PHI_Counter + 1 : PHI_Counter + 3, 1) =...
                                [OCPEC.p(i);...
                                OCPEC.w(i);...
                                OCPEC.symVar.s - OCPEC.p(i) * OCPEC.w(i)];
                            PHI_Counter = PHI_Counter + 3;
                        else
                            % box constraint variation inequality
                            OCPEC.PHI(PHI_Counter + 1 : PHI_Counter + 4, 1) =...
                                [OCPEC.p(i) - OCPEC.plant.l(i);...
                                OCPEC.plant.u(i) - OCPEC.p(i);...
                                OCPEC.symVar.s - (OCPEC.p(i) - OCPEC.plant.l(i)) * OCPEC.w(i);...
                                OCPEC.symVar.s + (OCPEC.plant.u(i) - OCPEC.p(i)) * OCPEC.w(i)];
                            PHI_Counter = PHI_Counter + 4;
                        end
                    end
            end
            
            %% initialize properties: dual variable
            % extend Dim struct
            OCPEC.Dim.sigma  = size(OCPEC.G, 1);            
            OCPEC.Dim.eta    = size(OCPEC.C, 1);
            OCPEC.Dim.lambda = OCPEC.Dim.x;
            OCPEC.Dim.gamma  = size(OCPEC.PHI, 1);  
            
            % initialize dual variable
            OCPEC.sigma  = sym('sigma', [OCPEC.Dim.sigma, 1]);            
            OCPEC.eta    = sym('eta', [OCPEC.Dim.eta, 1]);
            OCPEC.lambda = sym('lambda', [OCPEC.Dim.lambda, 1]);
            OCPEC.gamma  = sym('gamma', [OCPEC.Dim.gamma, 1]);             
           
            %% initialize properties: Feasibility Restoration Phase
            OCPEC.FRP.ZRef = zeros(OCPEC.Dim.Z, OCPEC.nStages);
            OCPEC.FRP.ZWeight = zeros(OCPEC.Dim.Z, OCPEC.nStages);
            
            OCPEC.symVar.ZRef = sym('ZRef', [OCPEC.Dim.Z, 1]);
            OCPEC.symVar.ZWeight = sym('ZWeight', [OCPEC.Dim.Z, 1]);
            Z = [OCPEC.tau; OCPEC.x; OCPEC.p; OCPEC.w];
            OCPEC.LFRP = 0.5 *(Z - OCPEC.symVar.ZRef)' * diag(OCPEC.symVar.ZWeight) * (Z - OCPEC.symVar.ZRef);
            
        end
        
    end
    
    %% Other Methods for OCPEC_Perturbed
    methods
        % set inequality and equality constraints by user defined symbolic function
        setInequalityConstraints(OCPEC, G)        
        
        setEqualityConstraints(OCPEC, C)  
        
        % code generation
        codeGen(OCPEC)
        
        % show OCPEC informulation
        showInfo(OCPEC)      
        
        % function and jacobian evaluation about cost and constraint
        L = computeCost_Function(OCPEC, Iterate, mode)
        
        Lvar = computeCost_Jacobian(OCPEC, Iterate, mode)
        
        [G, C, F] = computeConstraint_Function_G_C_F(OCPEC, Iterate)
        
        PHI = computeConstraint_Function_PHI(OCPEC, Iterate, s)
        
        [Gvar, Cvar, Fvar] = computeConstraint_Jacobian_G_C_F(OCPEC, Iterate)   
        
        PHIvar = computeConstraint_Jacobian_PHI(OCPEC, Iterate, s)   
            
    end
    
end

