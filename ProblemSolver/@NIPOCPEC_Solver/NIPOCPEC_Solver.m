classdef NIPOCPEC_Solver < handle
    %% NIPOCPEC_Solver: Superclass for applying the non interior point method to solve OCPEC
    %   We use FB function to project KKT conditions into a system of equation
    
    properties
        OCPEC; % obj, OCPEC problem
        Option; % struct, solver options 
    end
    
    %% Constructor Method for NIPOCPEC_Solver    
    methods
        function solver = NIPOCPEC_Solver(OCPEC)
            %NIPOCPEC_Solver
            %   Constructor of Class NIPOCPEC_Solver
            %
            % Syntax: 
            %           solver = NIPOCPEC_Solver(OCPEC)
            %
            % Argument: 
            %           OCPEC: obj, OCPEC problem
            %
            % Output: 
            %           solver: a object of Class NIPOCPEC_Solver   
            %
            %% construct an object
            solver.OCPEC = OCPEC;
            solver.Option = createSolverOption(solver); 
            
        end
        
    end
    
    %% Ohter Methods for NIPOCPEC_Solver    
    methods
        %% 
        Option = createSolverOption(solver)  
        
        checkOption(solver) % To Be Done
        
        showInfo(solver) 
        
        generateInitialGuess(solver)

        [solution, Info] = solveOCPEC(solver, IterateInit)
        
        showResult(solver, Info)
        
        %% Methods in solveOCPEC method
        % function and jacobian evaluation
        codeGen(solver)
        
        FunEval = FunctionEvaluation(solver, Iterate, s, z, mode)
        
        Hessian = computeHessian(solver, Iterate, s, mode)
        
        KKT_Residual = computeKKT_Residual(solver, Iterate, FunEval)
        
        KKT_Error = computeKKT_Error(solver, Iterate, FunEval, KKT_Residual)
        
        KKT_Matrix = computeKKT_Matrix(solver, FunEval)
        
        % FB
        [PSI, PSIa, PSIb] = computeFB_Function_Jacobian(solver, a, b, z)
        
        value = computeFB_minusInvPSIbPSI(solver, PSIb, PSI)
        
        % evaluate search direction
        [dY, Info] = SearchDirection_Riccati(solver, KKT_Residual, KKT_Matrix)       
        
        % Line Search     
        [Iterate_LS, Info] = LineSearch_Merit(solver, Iterate, FunEval, KKT_Residual, KKT_Matrix, beta, s, z, dY_k, mode)
        
        % Feasibility Restoration Phase
        [Iterate_FRP, Info] = FeasibilityRestorationPhase_MinDeviation(solver, Iterate_Ref, FunEval_Ref, s, z)
        
        % compute perturbed parameter s and z
        [s_k, z_k, FunEval_k] = computePerturedParam(solver, Iterate_k, FunEval_k, s, z)
 
        % examine solution
        Info = solutionExaminer(solver, solution, Record)
        
    end
    
end

