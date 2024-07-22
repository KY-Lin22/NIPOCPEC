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
        function self = NIP_C_Solver(OCPEC, NLP, Option)
            %NIP_C_Solver Construct an instance of this class
            %   Detailed explanation goes here         
            disp('creating solver...')
            % properties: OCPEC, NLP, Option, FuncObj
            self.OCPEC = OCPEC;
            self.NLP = NLP;     
            self.Option = Option;
            self.FuncObj = self.createFuncObj();
            disp('Done!')
        end
    end
    
    %% Other method
    methods
        % initialize properties              
        FuncObj = createFuncObj(self) 

        % main function of NIP_C
        [z_Opt, Info] = solveNLP(self, z_Init)

        % preprocess initial guess z_Init
        z_Init = preprocessInitialGuess(self, z_Init)

        % create parameter sequence
        [P, l_Max] = createParameterSequence(self)
     
        % main function of non-interior-point method
        [Y_Opt, Info] = NonInteriorPointMethod(self, z_Init, p_Init)

        % merit line search
        [Y_k, Info] = LineSearchMerit(self, beta, Y, p, dY)

        % Euler-Newton (predictor-corrector) continuation method
        [Y_l, Info] = EulerNewtonContinuationMethod(self, Y, p_l, p);

        % evaluate KKT error
        [KKT_error_primal, KKT_error_dual, KKT_error_complementary, KKT_error_total] = evaluateKKT_Error(self, Y, p)

        % evaluate natural residual
        natRes = evaluateNaturalResidual(self, z_Opt)
        
    end

    methods(Static)
        % create solver option
        Option = createSolverOption() 
    end
end

