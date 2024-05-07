classdef NIP_C_Solver < handle
    %Implementation of Non-Interior-Point Continuation (NIP_C) method which solves the following NLP problem:
    %  min  J(z),
    %  s.t. h(z) = 0,
    %       c(z) >= 0,
    %       g(z, s) >= 0
    % where: z is the variable, 
    %        s is the parameter, 
    %        J is the cost, and h, c, g are the constraint
    %    
    properties
        OCPEC % struct, optimal control problem with equalibrium constraints
        NLP % struct, nonlinear programming problem (discretized OCPEC)
        Option % struct, NIP_C solver option
        FuncObj % struct, function object
    end
    %% Constructor method
    methods
        function self = NIP_C_Solver(OCPEC, NLP)
            %NIP_C_Solver Construct an instance of this class
            %   Detailed explanation goes here
            
            disp('creating solver...')
            % properties: OCPEC and NLP
            self.OCPEC = OCPEC;
            self.NLP = NLP;     
            
            % properties: solver option
            self.Option = self.createSolverOption();
            
            % properties: function object
            self.FuncObj = self.createFuncObj();

            disp('Done!')
        end
    end
    
    %% Other method
    methods
        % initialize properties
        Option = createSolverOption(self)
        
        FuncObj = createFuncObj(self) 

        % main function of NIP_C
        [z_Opt, Info] = solveNLP(self, z_Init, p_Init, p_End) 

        % preprocess initial guess z_Init
        z_Init = preprocessInitialGuess(self, z_Init)

        % evaluate KKT matrix and sensitivity matrix
        KKT_Matrix = evaluateKKT_Matrix(self, h_grad, c_grad, g_grad, LAG_hessian,...
            PSI_c_grad_dual, PSI_c_grad_ineq, PSI_g_grad_dual, PSI_g_grad_ineq);
        
        Sensitivity_Matrix = evaluateSensitivity_Matrix(self,...
            PSI_c_grad_sigma, PSI_g_grad_ineq, PSI_g_grad_sigma)

        % evaluate KKT error
        [KKT_error_primal, KKT_error_dual, KKT_error_dual_scaled,...
            KKT_error_complementary, KKT_error_complementary_scaled, KKT_error_total] = ...
            evaluateKKT_error(self, gamma_h, gamma_c, gamma_g, h, c, g, LAG_grad_z)

        % evaluate natural residual
        natRes = evaluateNaturalResidual(self, z_Opt)       
        
        %% stage 1: Non-Interior-Point Method
        % main function of non-interior-point method
        [z_Opt, Info] = non_interior_point_method(self, z_Init, p)

        % merit line search
        [gamma_h_k, gamma_c_k, gamma_g_k, z_k, Info] = LineSearch_Merit(self, beta,...
            gamma_h, gamma_c, gamma_g, z, dgamma_h, dgamma_c, dgamma_g, dz,...
            s, sigma, J, h, PSI_c, PSI_g, J_grad)
        
        %% stage 2: Newton-Euler Continuation Method
        % Euler predict step
        [gamma_h_Eu_j, gamma_c_Eu_j, gamma_g_Eu_j, z_Eu_j, Info] = Euler_predict_step(self,...
            gamma_h, gamma_c, gamma_g, z, p, p_j)
 
        % Newton correct step
        [gamma_h_j, gamma_c_j, gamma_g_j, z_j, Info] = Newton_correct_step(self,...
            gamma_h_Eu_j, gamma_c_Eu_j, gamma_g_Eu_j, z_Eu_j, p_j)
    end
end

