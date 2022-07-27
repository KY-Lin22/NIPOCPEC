function [Iterate_FRP, Info] = FeasibilityRestorationPhase_MinDeviation(solver, Iterate_Ref, FunEval_Ref, s, z)
%FeasibilityRestorationPhase_MinDeviation
%   Detailed explanation goes here

TimeStart = tic;
switch solver.Option.printLevel
    case 0
    case 1
        disp(' ')
        disp('*---- Entering Feasibility Restoration Phase...----*')
    case 2
        disp(' ')
        disp('*---- Entering Feasibility Restoration Phase...----*')
end

% compute ZRef and ZWeight
ZRef = [Iterate_Ref.tau; Iterate_Ref.x; Iterate_Ref.p; Iterate_Ref.w];
ZWeight = zeros(solver.OCPEC.Dim.Z, solver.OCPEC.nStages);
nu_Z = solver.Option.FRP.nu_Z;
for n = 1 : solver.OCPEC.nStages
    ZWeight(:, n) = nu_Z * min([ones(solver.OCPEC.Dim.Z, 1), 1./abs(ZRef(:, n))], [], 2);    
end
solver.OCPEC.FRP.ZRef = ZRef;
solver.OCPEC.FRP.ZWeight = ZWeight;

%%
OCPEC = solver.OCPEC;
Option = solver.Option;
nStages = OCPEC.nStages;
Dim = OCPEC.Dim;
VI_mode = OCPEC.VI_mode;
maxIterNum = Option.FRP.maxIterNum;
nu_M = Option.FRP.nu_M;
betaInit = Option.FRP.betaInit;
employLSMN = Option.FRP.employLeastSquareMinNorm;
% create IterateInit (reusing Iterate_Ref except the dual variables of equality-type constraint)
IterateInit = struct('tau', Iterate_Ref.tau, 'x', Iterate_Ref.x, 'p', Iterate_Ref.p, 'w', Iterate_Ref.w,...
    'sigma', Iterate_Ref.sigma, 'eta', zeros(size(Iterate_Ref.eta)), 'lambda', zeros(size(Iterate_Ref.lambda)), 'gamma', []);
if (strcmp(VI_mode,'Reg_NCPs')) || (strcmp(VI_mode, 'Reg_Scholtes'))
    IterateInit.gamma = Iterate_Ref.gamma;
elseif strcmp(VI_mode, 'SmoothingEquation')
    IterateInit.gamma = zeros(size(Iterate_Ref.gamma));
end

%% solving OCPEC
% initialize FRP iteration routine (x: previous iterate; x_j: current iterate)
Iterate = IterateInit;
beta = betaInit;

FailureFlag.LineSearch = false;

% FRP iteration routine 
for j = 1 : maxIterNum + 1
    %% step 1: Function and Jacobian Evaluation of Previous Iterate (totalCost and Feasibility)    
    % function
    if j == 1
        % only the first one needs(reusing FunEval_Ref except the cost function), others can reuse the results from line search       
        FunEval = struct('L', [],...
            'G', FunEval_Ref.G, 'C', FunEval_Ref.C, 'F', FunEval_Ref.F, 'PHI',FunEval_Ref.PHI,...
            'PSIg', FunEval_Ref.PSIg, 'PSIgSigma', FunEval_Ref.PSIgSigma, 'PSIgG', FunEval_Ref.PSIgG,...
            'PSIphi', FunEval_Ref.PSIphi, 'PSIphiGamma', FunEval_Ref.PSIphiGamma, 'PSIphiPHI', FunEval_Ref.PSIphiPHI,...
            'Lvar', [], 'Gvar', [], 'Cvar', [], 'Fvar', [], 'PHIvar',[],...
            'Hessian', []);
        FunEval.L = OCPEC.computeCost_Function(Iterate, 'FRP');
    end    
    % totalCost and Feasibility
    totalCost = sum(FunEval.L);
    if (strcmp(VI_mode,'Reg_NCPs')) || (strcmp(VI_mode, 'Reg_Scholtes'))
        Feasibility = norm(reshape([FunEval.PSIg; FunEval.C; FunEval.F; FunEval.PSIphi], [], 1), 1);
    elseif strcmp(VI_mode, 'SmoothingEquation')
        Feasibility = norm(reshape([FunEval.PSIg; FunEval.C; FunEval.F; FunEval.PHI], [], 1), 1);
    end  
    
    %% Record and Print Information of Previous Iterate   
    % record
    if j == 1
        totalCostInit = totalCost;
        FeasibilityInit = Feasibility;
    end      
    % print
    switch Option.printLevel
        case 0
            prevIterMsg = [];
        case 1
            prevIterMsg = [];
        case 2
            prevIterMsg = ['Iter: ', num2str(j - 1), '; ',...
                'Cost: ', num2str(totalCost, '%10.3e'), '; ',...
                'Feasibility: ', num2str(Feasibility, '%10.3e'), '; ',...
                'beta: ' num2str(beta, '%10.3e')];
    end
    disp(prevIterMsg)    
    
    %% step 2: Checking Termination
    terminalCond.Feasibility = (Feasibility <= nu_M * FeasibilityInit);    
    terminalCond.maxIterNum = (j == (maxIterNum + 1));
    terminalCond.LS = FailureFlag.LineSearch;  
    
    if terminalCond.Feasibility
        % FRP finds a less infeasibility solution
        exitFlag = true;
        solutionType = 1;
        FRP_FailureFlag = false;
    elseif terminalCond.maxIterNum || terminalCond.LS
        % FRP fails to find the a less infeasibility solution
        exitFlag = true;
        solutionType = 2;
        FRP_FailureFlag = true;
    else
        exitFlag = false;        
    end        
    
    %% step 3: Checking Exit Flag
    if exitFlag
        % return solution
        switch solutionType
            case 1
                % return a less infeasibility solution (reusing Iterate except the dual variables of equality-type constraint)
                Iterate_FRP = struct('tau', Iterate.tau, 'x', Iterate.x, 'p', Iterate.p, 'w', Iterate.w,...
                    'sigma', Iterate.sigma, 'eta', zeros(Dim.eta, nStages), 'lambda', zeros(Dim.lambda, nStages), 'gamma', []);                
                if (strcmp(VI_mode, 'Reg_NCPs')) || (strcmp(VI_mode, 'Reg_Scholtes'))
                    Iterate_FRP.gamma  = Iterate.gamma;
                elseif (strcmp(VI_mode, 'SmoothingEquation'))
                    Iterate_FRP.gamma  = zeros(Dim.gamma, nStages);
                end
                % compute the dual variables of equality-type constraint using lsqminnorm (Optional)          
                if employLSMN
                    [Iterate_FRP.eta, Iterate_FRP.lambda, Iterate_FRP.gamma] = computeDualVar_lsqminnorm(solver, Iterate_FRP, s);
                end               
                
                % create FunEval_FRP (reusing FunEval except cost function)
                FunEval_FRP = struct('L', [],...
                    'G', FunEval.G, 'C', FunEval.C, 'F', FunEval.F, 'PHI',FunEval.PHI,...
                    'PSIg', FunEval.PSIg, 'PSIgSigma', FunEval.PSIgSigma, 'PSIgG', FunEval.PSIgG,...
                    'PSIphi', FunEval.PSIphi, 'PSIphiGamma', FunEval.PSIphiGamma, 'PSIphiPHI', FunEval.PSIphiPHI,...
                    'Lvar', [], 'Gvar', [], 'Cvar', [], 'Fvar', [], 'PHIvar',[],...
                    'Hessian', []);
                FunEval_FRP.L = OCPEC.computeCost_Function(Iterate_FRP, 'Regular');                  
                
                % termination message
                terminationMsg_1 = ['FRP returns a less infeasibility iterate after ', num2str(j - 1), ' iteration; '];
                terminationMsg_2 = ['Feasibility: ', num2str(FeasibilityInit), '(Init)-->', num2str(Feasibility), '(End); ',...
                    'DeviationCost: ' num2str(totalCostInit), '(Init)-->', num2str(totalCost), '(End); '];                
            case 2
                % return Iterate_Ref and FunEval_Ref
                Iterate_FRP = Iterate_Ref;
                FunEval_FRP = FunEval_Ref;                
                
                % termination message
                terminationMsg_1 = 'FRP fails to return a less infeasibility iterate';
                terminationMsg_2 = ['failure condition are: ',...
                    'maxIterNum', '(', mat2str(terminalCond.maxIterNum), '); ', 'LineSearch', '(', mat2str(terminalCond.LS), '):'];        
        end
        % record and print information about solution
        FRP_TimeElapsed = toc(TimeStart);
        Info.FunEval = FunEval_FRP;
        Info.Time = FRP_TimeElapsed;
        Info.FailureFlag = FRP_FailureFlag;
        
        switch Option.printLevel
            case 0
            case 1
                disp(terminationMsg_1)
                disp(terminationMsg_2)
            case 2
                disp(terminationMsg_1)
                disp(terminationMsg_2)
        end

        break
    end
    
    %% step 4: Function and Jacobian Evaluation of Previous Iterate (KKT Residual and Matrix)
    % Jacobian and KKT residual    
    if j == 1
        % reusing FunEval_Ref except the cost function
        FunEval.Gvar = FunEval_Ref.Gvar;
        FunEval.Cvar = FunEval_Ref.Cvar;
        FunEval.Fvar = FunEval_Ref.Fvar;
        FunEval.PHIvar = FunEval_Ref.PHIvar;
    else
        [FunEval.Gvar, FunEval.Cvar, FunEval.Fvar] = OCPEC.computeConstraint_Jacobian_G_C_F(Iterate);
        FunEval.PHIvar = OCPEC.computeConstraint_Jacobian_PHI(Iterate, s);
    end
    FunEval.Lvar = OCPEC.computeCost_Jacobian(Iterate, 'FRP');
    
    KKT_Residual = solver.computeKKT_Residual(Iterate, FunEval);
    
    % Hessian and KKT matrix 
    FunEval.Hessian = solver.computeHessian(Iterate, s, 'FRP');
    KKT_Matrix = solver.computeKKT_Matrix(FunEval);
    
    %% step 5: Search Direction Evaluation
    [dY_k, ~] = solver.SearchDirection_Riccati(KKT_Residual, KKT_Matrix);    
    
    %% step 6: Merit Line Search
    [Iterate_LS, Info_LS] = solver.LineSearch_Merit(Iterate, FunEval, KKT_Residual, KKT_Matrix, beta, s, z, dY_k, 'FRP');
    % check failure flag and determine new iterate and its function evaluation
    if Info_LS.FailureFlag
        FailureFlag.LineSearch = true;
        % use previous iterate       
        Iterate_k  = Iterate;
        beta_k     = beta;
        FunEval_k  = FunEval;        
    else
        FailureFlag.LineSearch = false;
        % use MLS iterate     
        Iterate_k  = Iterate_LS;
        beta_k     = Info_LS.beta;
        FunEval_k  = Info_LS.FunEval;        
    end    
    % prepare for next iteration
    Iterate = Iterate_k;
    beta = beta_k;
    FunEval = FunEval_k;
end

switch solver.Option.printLevel
    case 0
    case 1
        disp('*--------------------------------------------------*')
    case 2
        disp('*--------------------------------------------------*')
end

end

%%
function [eta, lambda, gamma] = computeDualVar_lsqminnorm(solver, Iterate, s)
OCPEC = solver.OCPEC;
VI_mode = OCPEC.VI_mode;
Dim = OCPEC.Dim;
nStages = OCPEC.nStages;
LAMBDA_threshold = 1000;

% Jacobian evaluation
Lvar = OCPEC.computeCost_Jacobian(Iterate, 'Regular');
[Gvar, Cvar, Fvar] = OCPEC.computeConstraint_Jacobian_G_C_F(Iterate);
PHIvar = OCPEC.computeConstraint_Jacobian_PHI(Iterate, s);

%% compute dual variables by solving an overdetermined linear equation system in a backward recursion manner
eta = zeros(Dim.eta, nStages);
lambda = zeros(Dim.lambda, nStages);
gamma = zeros(Dim.gamma, nStages);
% initialize lambdaNext_n
lambdaNext_n = zeros(Dim.lambda, 1);
%
for n = nStages: -1 : 1
    % load
    Ltau_n = Lvar.Ltau(:, 1 + (n - 1) * Dim.tau : n * Dim.tau);
    Lx_n   = Lvar.Lx(:, 1 + (n - 1) * Dim.x : n * Dim.x);
    Lp_n   = Lvar.Lp(:, 1 + (n - 1) * Dim.p : n * Dim.p);
    Lw_n   = Lvar.Lw(:, 1 + (n - 1) * Dim.w : n * Dim.w);
    
    Gtau_n = Gvar.Gtau(:, 1 + (n - 1) * Dim.tau : n * Dim.tau);
    Gx_n   = Gvar.Gx(:, 1 + (n - 1) * Dim.x : n * Dim.x);
    Gp_n   = Gvar.Gp(:, 1 + (n - 1) * Dim.p : n * Dim.p);    
    
    Ctau_n = Cvar.Ctau(:, 1 + (n - 1) * Dim.tau : n * Dim.tau);
    Cx_n   = Cvar.Cx(:, 1 + (n - 1) * Dim.x : n * Dim.x);
    Cp_n   = Cvar.Cp(:, 1 + (n - 1) * Dim.p : n * Dim.p);
    Cw_n   = Cvar.Cw(:, 1 + (n - 1) * Dim.w : n * Dim.w);
    Cvar_n = [Ctau_n, Cx_n, Cp_n, Cw_n];

    Ftau_n = Fvar.Ftau(:, 1 + (n - 1) * Dim.tau : n * Dim.tau);
    Fx_n   = Fvar.Fx(:, 1 + (n - 1) * Dim.x : n * Dim.x);
    Fp_n   = Fvar.Fp(:, 1 + (n - 1) * Dim.p : n * Dim.p);
    Fvar_n = [Ftau_n, Fx_n, Fp_n, zeros(Dim.lambda, Dim.w)];
    
    PHIp_n = PHIvar.PHIp(:, 1 + (n - 1) * Dim.p : n * Dim.p);
    PHIw_n = PHIvar.PHIw(:, 1 + (n - 1) * Dim.w : n * Dim.w);
    
    % solve overdetermined linear equation system
    sigma_n = Iterate.sigma(:, n);   
    if (strcmp(VI_mode,'Reg_NCPs')) || (strcmp(VI_mode, 'Reg_Scholtes'))
        gamma_n = Iterate.gamma(:, n);
        % overdetermined formulate linear system
        A_n = [Cvar_n', Fvar_n'];
        b_n = [Ltau_n' - Gtau_n' * sigma_n;...
               Lx_n'   - Gx_n'   * sigma_n                      + lambdaNext_n;...
               Lp_n'   - Gp_n'   * sigma_n  - PHIp_n' * gamma_n;...
               Lw_n'                        - PHIw_n' * gamma_n];
        % solve linear system and extract eta, lambda
        EtaLambda_n = lsqminnorm(A_n,-b_n);
        if norm(EtaLambda_n, Inf) > LAMBDA_threshold
            EtaLambda_n = zeros(Dim.eta + Dim.lambda, 1);
        end
        eta_n    = EtaLambda_n(1 : Dim.eta, :);
        lambda_n = EtaLambda_n(Dim.eta + 1 : end, :);
    elseif strcmp(VI_mode, 'SmoothingEquation')
        % overdetermined formulate linear system
        PHIvar_n = [zeros(Dim.gamma, Dim.tau + Dim.x), PHIp_n, PHIw_n];
        A_n = [Cvar_n', Fvar_n', PHIvar_n'];
        b_n = [Ltau_n' - Gtau_n' * sigma_n;...
               Lx_n'   - Gx_n'   * sigma_n  + lambdaNext_n;...
               Lp_n'   - Gp_n'   * sigma_n;...
               Lw_n'    ];
        % solve linear system and extract eta, lambda, gamma
        EtaLambdaGamma_n = lsqminnorm(A_n,-b_n);
        if norm(EtaLambdaGamma_n, Inf) > LAMBDA_threshold
            EtaLambdaGamma_n = zeros(Dim.eta + Dim.lambda + Dim.gamma, 1);
        end
        eta_n    = EtaLambdaGamma_n(1 : Dim.eta, :);
        lambda_n = EtaLambdaGamma_n(Dim.eta + 1 : Dim.eta + Dim.lambda, :);
        gamma_n  = EtaLambdaGamma_n(Dim.eta + Dim.lambda + 1 : end, :);       
    end
    
    % record and prepare for next iteration
    eta(:, n) = eta_n;
    lambda(:, n) = lambda_n;
    gamma(:, n) = gamma_n;
    lambdaNext_n = lambda_n; 
end

end
