function [solution, Info] = solveOCPEC(solver, IterateInit)
%solveOCPEC
%   Detailed explanation goes here

OCPEC = solver.OCPEC;
Option = solver.Option;
maxIterNum = Option.maxIterNum;
Tol = Option.Tolerance;
betaInit = Option.LineSearch.betaInit;
employFRP = Option.employFeasibilityRestorationPhase;

sInit = Option.sInit;
zInit = Option.zInit;
sEnd = Option.sEnd;
zEnd = Option.zEnd;

% create record
Record = struct('s', zeros(maxIterNum + 1, 1), 'z', zeros(maxIterNum + 1, 1), 'IterateType', [],...
    'Merit', zeros(maxIterNum + 1, 1), 'beta', zeros(maxIterNum + 1, 1), 'stepSize', zeros(maxIterNum + 1, 1), ...
    'totalCost', zeros(maxIterNum + 1, 1), 'KKT_Error', [], 'Time', [],...
    'iterNum', 0, 'terminalStatus',0, 'terminalCond', []); 
Record.IterateType = cell(maxIterNum + 1, 1);
Record.KKT_Error = struct('Total', zeros(maxIterNum + 1, 1), 'Feasibility', zeros(maxIterNum + 1, 1), 'Stationarity', zeros(maxIterNum + 1, 1));   
Record.Time = struct('FunEval', 0, 'KKT', 0, 'SearchDirection', 0, 'LineSearch', 0, 'FRP', 0, 'else', 0, 'total', 0);

%% solving OCPEC
disp('******************************************************************');
disp('Computing the Optimal Solution for OCPEC...');

% initialize regular iteration routine (x: previous iterate; x_k: current iterate)
s = sInit;
z = zInit;
Iterate = IterateInit;
IterateType = 'Init';
Merit = 0;
beta = betaInit;
stepSize = 0;

TimeElasped_FunEval = 0;
TimeElasped_KKT = 0;
TimeElasped_SearchDirection = 0;
TimeElasped_LineSearch = 0;
TimeElasped_FRP = 0;
TimeElasped_Total = 0;

FailureFlag.LineSearch = false;
FailureFlag.FRP = false;

% regular iteration routine
for k = 1 : maxIterNum + 1
    totalTimeStart = tic;
    %% step 1: Function and Jacobian Evaluation of Previous Iterate (KKT Residual)     
    FunEval_Jacobian_TimeStart = tic;
    
    % function
    if k == 1
        % only the first one needs, others can reuse the results from line search or FRP
        FunEval = solver.FunctionEvaluation(Iterate, s, z, 'Regular');
    end    
    % Jacobian
    FunEval.Lvar = OCPEC.computeCostFunJacobian(Iterate, 'Regular');
    [FunEval.Gvar, FunEval.Cvar, FunEval.Fvar] = OCPEC.computeConstraintFunJacobian_G_C_F(Iterate);
    FunEval.PHIvar = OCPEC.computeConstraintFunJacobian_PHI(Iterate, s);   
    
    TimeElasped_FunEval_Jacobian = toc(FunEval_Jacobian_TimeStart);
    
    % totalCost, KKT residual and error
    KKT_Residual_TimeStart = tic;
    
    totalCost = sum(FunEval.L);
    KKT_Residual = solver.computeKKT_Residual(Iterate, FunEval);
    KKT_Error = solver.computeKKT_Error(Iterate, KKT_Residual);
    
    TimeElasped_KKT_Residual = toc(KKT_Residual_TimeStart);
    
    %% Record and Print Information of Previous Iterate
    Record.s(k) = s;
    Record.z(k) = z;
    Record.IterateType{k} = IterateType;
    
    Record.Merit(k) = Merit;
    Record.beta(k) = beta;
    Record.stepSize(k) = stepSize;    
    
    Record.totalCost(k)              = totalCost;
    Record.KKT_Error.Total(k)        = KKT_Error.Total;
    Record.KKT_Error.Feasibility(k)  = KKT_Error.Feasibility;
    Record.KKT_Error.Stationarity(k) = KKT_Error.Stationarity;  
    
    TimeElasped_Else = TimeElasped_Total - TimeElasped_FunEval - TimeElasped_KKT...
        - TimeElasped_SearchDirection - TimeElasped_LineSearch - TimeElasped_FRP;    
    Record.Time.FunEval         = Record.Time.FunEval         + TimeElasped_FunEval;
    Record.Time.KKT             = Record.Time.KKT             + TimeElasped_KKT;
    Record.Time.SearchDirection = Record.Time.SearchDirection + TimeElasped_SearchDirection;
    Record.Time.LineSearch      = Record.Time.LineSearch      + TimeElasped_LineSearch;
    Record.Time.FRP             = Record.Time.FRP             + TimeElasped_FRP;    
    Record.Time.else            = Record.Time.else            + TimeElasped_Else;
    Record.Time.total           = Record.Time.total           + TimeElasped_Total;
    
    switch Option.printLevel
        case 0
            prevIterMsg = [];
        case 1
            prevIterMsg = [];
        case 2
            prevIterMsg = ['Iter: ', num2str(k - 1), '(', IterateType, ')', '; ',...
                's: ' num2str(s,'%10.2e'), '; ',...
                'z: ', num2str(z,'%10.2e'), '; ',...
                'Cost: ', num2str(totalCost,'%10.2e'), '; ',...
                'KKT: ',num2str(KKT_Error.Feasibility,'%10.2e'), '(F) ',num2str(KKT_Error.Stationarity,'%10.2e'),'(S); ',...
                'Merit: ', num2str(Merit,'%10.2e'), '; ',...
                'beta: ' num2str(beta,'%10.2e'), '; ',...
                'StepSize: ', num2str(stepSize,'%10.2e'), '; ', ...
                'Time: ', num2str(1000 * TimeElasped_Total,'%10.2e'), ' ms'];           
    end
    disp(prevIterMsg);
    
    %% step 2: Checking Termination
    terminalCond.sz = ((s == sEnd) && (z == zEnd));
    terminalCond.KKT_T = (KKT_Error.Total <= Tol.KKT_Error_Total);
    terminalCond.KKT_F = (KKT_Error.Feasibility <= Tol.KKT_Error_Feasibility); 
    terminalCond.KKT_S = (KKT_Error.Stationarity <= Tol.KKT_Error_Stationarity);     
    terminalCond.maxIterNum = (k == (maxIterNum + 1));
    terminalCond.LSwithoutFRP = (FailureFlag.LineSearch && ~employFRP);
    terminalCond.LSwithFRP = (FailureFlag.LineSearch && FailureFlag.FRP);    
    
    if  terminalCond.sz && (terminalCond.KKT_T || terminalCond.KKT_F || terminalCond.KKT_S)
        % solver finds the optimal solution
        exitFlag = true;
        terminalStatus = 1;
    elseif terminalCond.maxIterNum || terminalCond.LSwithoutFRP || terminalCond.LSwithFRP
        % solver fails to find the optimal solution
        exitFlag = true;
        terminalStatus = 0;
    else
        exitFlag = false;
    end
    
    %% step 3: Checking Exit Flag
    if exitFlag
        % return solution (i.e.,iterate)
        Record.iterNum = k - 1; 
        Record.terminalStatus = terminalStatus;     
        Record.terminalCond = terminalCond;
        solution = Iterate; 
        Info = solver.solutionExaminer(solution, Record);% haven't done
        break
    end
       
    %% step 4: Function and Jacobian Evaluation of Previous Iterate (KKT Matrix)
    % Hessian and KKT matrix 
    FunEval_Hessian_TimeStart = tic;
    
    FunEval.Hessian = solver.computeHessian(Iterate, s, 'Regular');
    
    TimeElasped_FunEval_Hessian = toc(FunEval_Hessian_TimeStart);
    TimeElasped_FunEval = TimeElasped_FunEval_Jacobian + TimeElasped_FunEval_Hessian;
    
    % KKT matrix
    KKT_Matrix_TimeStart = tic;
    
    KKT_Matrix = solver.computeKKT_Matrix(FunEval);
    
    TimeElasped_KKT_Matrix = toc(KKT_Matrix_TimeStart);
    TimeElasped_KKT = TimeElasped_KKT_Residual + TimeElasped_KKT_Matrix;  
    
    %% step 5: Search Direction Evaluation
    [dY_k, Info_SearchDirection] = solver.SearchDirection_Riccati(KKT_Residual, KKT_Matrix);
    TimeElasped_SearchDirection = Info_SearchDirection.Time;    
    
    %% step 6: Merit Line Search
    [Iterate_LS, Info_LS] = solver.LineSearch_Merit(Iterate, FunEval, KKT_Residual, KKT_Matrix, beta, s, z, dY_k, 'Regular');
    TimeElasped_LineSearch = Info_LS.Time;
    
    % check failure flag
    if Info_LS.FailureFlag
        FailureFlag.LineSearch = true;
        IterateType_k = 'Prev';       
    else
        FailureFlag.LineSearch = false;
        IterateType_k = Info_LS.IterateType;      
    end
    
    %% Step 7: Feasibility Restoration Phase
    if FailureFlag.LineSearch && employFRP
        [Iterate_FRP, Info_FRP] = solver.FeasibilityRestorationPhase_MinDeviation(Iterate, FunEval, s, z);
        TimeElasped_FRP = Info_FRP.Time;
        % check failure flag
        if Info_FRP.FailureFlag
            FailureFlag.FRP = true;
            IterateType_k = 'Prev';
        else
            FailureFlag.FRP = false;           
            FailureFlag.LineSearch = false; % reset FailureFlag due to the success of FRP
            IterateType_k = 'FRP';
        end        
    else       
        TimeElasped_FRP = 0;
    end    
    
    %% Step 8: Determine New Iterate and Its Function Evaluation
    if (strcmp(IterateType_k,'MLS')) || (strcmp(IterateType_k, 'SOC'))
        % use MLS or SOC iterate
        Iterate_k  = Iterate_LS;        
        Merit_k    = Info_LS.Merit;
        beta_k     = Info_LS.beta;
        stepSize_k = Info_LS.stepSize;
        FunEval_k  = Info_LS.FunEval;       
    elseif (strcmp(IterateType_k,'FRP'))
        % use FRP iterate
        Iterate_k  = Iterate_FRP;       
        Merit_k    = 0;
        beta_k     = betaInit;% reset penalty parameter
        stepSize_k = 1;
        FunEval_k  = Info_FRP.FunEval;          
    elseif (strcmp(IterateType_k,'Prev'))
        % use previous iterate
        Iterate_k  = Iterate;
        Merit_k    = Merit;
        beta_k     = beta;
        stepSize_k = stepSize;
        FunEval_k  = FunEval;
    else
        error('wrong iterate type')
    end
    
    %% Step 9: Update Perturbed Parameter and Prepare for Next Iteration
    % update perturbed parameter s, z and FunEval_k
    if (strcmp(IterateType_k,'MLS')) || (strcmp(IterateType_k, 'SOC'))
        [s_k, z_k, FunEval_k] = solver.computePerturedParam(Iterate_k, FunEval_k, s, z);
    elseif (strcmp(IterateType_k,'FRP')) || (strcmp(IterateType_k,'Prev'))
        s_k = s;
        z_k = z;        
    else
        error('wrong iterate type')
    end
    % reset penalty parameter
    if (s_k ~= s) || (z_k ~= z)        
        beta_k = betaInit;
    end    
    % prepare for next iteration
    s = s_k;
    z = z_k;
    Iterate = Iterate_k;
    IterateType = IterateType_k;
    Merit = Merit_k;
    beta = beta_k;
    stepSize = stepSize_k;
    
    FunEval = FunEval_k;
    
    TimeElasped_Total = toc(totalTimeStart);
end

disp('******************************************************************');

end

