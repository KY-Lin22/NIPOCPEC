function [Iterate_LS, Info] = LineSearch_Merit(solver, Iterate, FunEval, KKT_Residual, KKT_Matrix, beta, s, z, dY_k, mode)
%LineSearch_Merit
%   Detailed explanation goes here
TimeStart = tic;

OCPEC = solver.OCPEC;
Option = solver.Option;
Dim = OCPEC.Dim;
nStages = OCPEC.nStages;
employSOC = Option.employSecondOrderCorrection;

nu_D = Option.LineSearch.nu_D;
stepSize_Init = 1;
switch mode
    case 'Regular'
        rho = Option.LineSearch.rho;
        stepSize_Min = Option.LineSearch.stepSize_Min;
        stepSize_DecayRate = Option.LineSearch.stepSize_DecayRate;
    case 'FRP'
        rho = Option.FRP.rho;
        stepSize_Min = Option.FRP.stepSize_Min;
        stepSize_DecayRate = Option.FRP.stepSize_DecayRate;      
end

%% Some Evaluation Quantities of Previous Iterate
% cost and its directional derivative
totalCost = sum(FunEval.L);

dZ_k = dY_k(Dim.Node(4) + 1 : Dim.Node(8), :);
totalCostDD = 0;
for n = 1 : nStages 
    LZ_n = [FunEval.Lvar.Ltau(:, 1 + (n - 1) * Dim.tau : n * Dim.tau),...
        FunEval.Lvar.Lx(:, 1 + (n - 1) * Dim.x : n * Dim.x),...
        FunEval.Lvar.Lp(:, 1 + (n - 1) * Dim.p : n * Dim.p),...
        FunEval.Lvar.Lw(:, 1 + (n - 1) * Dim.w : n * Dim.w)];      
    totalCostDD = totalCostDD + LZ_n * dZ_k(:, n);
end

% constraint violation (L1 norm)
if (strcmp(OCPEC.VI_mode,'Reg_NCPs')) || (strcmp(OCPEC.VI_mode, 'Reg_Scholtes'))
    totalCstrVio_L1Norm = norm(reshape([FunEval.PSIg; FunEval.C; FunEval.F; FunEval.PSIphi], [], 1), 1); 
elseif strcmp(OCPEC.VI_mode, 'SmoothingEquation')
    totalCstrVio_L1Norm = norm(reshape([FunEval.PSIg; FunEval.C; FunEval.F; FunEval.PHI], [], 1), 1); 
end

% penalty parameter
beta_Trial = totalCostDD/((1 - rho) * totalCstrVio_L1Norm);
if beta >= beta_Trial
    beta_k = beta;
else
    beta_k = beta_Trial + 1;
end

% merit and its directional derivative 
merit = totalCost + beta_k * totalCstrVio_L1Norm;
meritDD = totalCostDD - beta_k * totalCstrVio_L1Norm;

%% Backtracking Line Search
hasFoundNewIterate = false;

while ~hasFoundNewIterate
    %% Step 1: estimate trail stepsize, iterate and merit 
    stepSize_trial = max([stepSize_Init, stepSize_Min]);
    
    Iterate_trial.sigma  = Iterate.sigma  + stepSize_trial * dY_k(              1 : Dim.Node(1), :);
    Iterate_trial.eta    = Iterate.eta    + stepSize_trial * dY_k(Dim.Node(1) + 1 : Dim.Node(2), :);
    Iterate_trial.lambda = Iterate.lambda + stepSize_trial * dY_k(Dim.Node(2) + 1 : Dim.Node(3), :);
    Iterate_trial.gamma  = Iterate.gamma  + stepSize_trial * dY_k(Dim.Node(3) + 1 : Dim.Node(4), :);
    Iterate_trial.tau    = Iterate.tau    + stepSize_trial * dY_k(Dim.Node(4) + 1 : Dim.Node(5), :);
    Iterate_trial.x      = Iterate.x      + stepSize_trial * dY_k(Dim.Node(5) + 1 : Dim.Node(6), :);
    Iterate_trial.p      = Iterate.p      + stepSize_trial * dY_k(Dim.Node(6) + 1 : Dim.Node(7), :);
    Iterate_trial.w      = Iterate.w      + stepSize_trial * dY_k(Dim.Node(7) + 1 : Dim.Node(8), :);
    
    FunEval_trial = solver.FunctionEvaluation(Iterate_trial, s, z, mode);
    if (strcmp(OCPEC.VI_mode,'Reg_NCPs')) || (strcmp(OCPEC.VI_mode, 'Reg_Scholtes'))
        totalCstrVio_L1Norm_trial = norm(reshape([FunEval_trial.PSIg; FunEval_trial.C; FunEval_trial.F; FunEval_trial.PSIphi], [], 1), 1);
    elseif strcmp(OCPEC.VI_mode,  'SmoothingEquation')
        totalCstrVio_L1Norm_trial = norm(reshape([FunEval_trial.PSIg; FunEval_trial.C; FunEval_trial.F; FunEval_trial.PHI], [], 1), 1);
    end
    merit_trial = sum(FunEval_trial.L) + beta_k * totalCstrVio_L1Norm_trial;
    
    %% Step 2: checking sufficient decrease condition
    if merit_trial <= merit + stepSize_trial * nu_D * meritDD
        % return merit line search iterate
        hasFoundNewIterate = true;
        LineSearchFailureFlag = false;
        IterateType = 'MLS';     
        
    elseif (strcmp(mode, 'Regular')) && (employSOC) && (stepSize_trial == 1) && (sum(FunEval_trial.L) <= totalCost) 
        % estimate second order correction (SOC) iterate and merit
        Iterate_trial = SecondOrderCorrection(solver, Iterate, FunEval, KKT_Residual, KKT_Matrix, FunEval_trial);
        FunEval_trial= solver.FunctionEvaluation(Iterate_trial, s, z, mode);
        if (strcmp(OCPEC.VI_mode,'Reg_NCPs')) || (strcmp(OCPEC.VI_mode, 'Reg_Scholtes'))
            totalCstrVio_L1Norm_trial = norm(reshape([FunEval_trial.PSIg; FunEval_trial.C; FunEval_trial.F; FunEval_trial.PSIphi], [], 1), 1);
        elseif strcmp(OCPEC.VI_mode,  'SmoothingEquation')
            totalCstrVio_L1Norm_trial = norm(reshape([FunEval_trial.PSIg; FunEval_trial.C; FunEval_trial.F; FunEval_trial.PHI], [], 1), 1);
        end
        merit_trial = sum(FunEval_trial.L) + beta_k * totalCstrVio_L1Norm_trial;
        
        if merit_trial <= merit + stepSize_trial * nu_D * meritDD
            % return SOC iterate
            hasFoundNewIterate = true;
            LineSearchFailureFlag = false;
            IterateType = 'SOC';            
        else
            % discard this SOC step and resume backtracking line search with smaller stepsize
            stepSize_Init = stepSize_DecayRate * stepSize_Init; 
        end
        
    else
        % need to estimate a smaller stepsize
        stepSize_Init = stepSize_DecayRate * stepSize_Init;
        
    end  
    
    %% Step 3: checking min stepsize
    if (stepSize_trial == stepSize_Min)&&(~hasFoundNewIterate)
        LineSearchFailureFlag = true;
        IterateType = 'MLS'; % return the iterate_trail obtained by stepSize_Min          
        break
    end     
    
end

%% Record Information
TimeElapsed = toc(TimeStart);

Iterate_LS = Iterate_trial;

Info.IterateType = IterateType;
Info.Merit = merit_trial;
Info.beta = beta_k;
Info.stepSize = stepSize_trial;
Info.FunEval = FunEval_trial;

Info.FailureFlag = LineSearchFailureFlag;
Info.Time = TimeElapsed;

end

%%
function Iterate_SOC = SecondOrderCorrection(solver, Iterate, FunEval, KKT_Residual, KKT_Matrix, FunEval_full)
%SecondOrderCorrection
%   Detailed explanation goes here

OCPEC = solver.OCPEC;
Dim = OCPEC.Dim;

% SOC KKT residual
KKT_Residual_SOC = struct('G_Fsb', [], 'C_Fsb', [], 'F_Fsb', [], 'PHI_Fsb', [],...
    'HtauT', KKT_Residual.HtauT, 'HxTlambdaNext', KKT_Residual.HxTlambdaNext, 'HpT', KKT_Residual.HpT, 'HwT', KKT_Residual.HwT);

G_Fsb_Correct = solver.computeFB_minusInvPSIbPSI(FunEval.PSIgG, FunEval_full.PSIg);
KKT_Residual_SOC.G_Fsb = KKT_Residual.G_Fsb + G_Fsb_Correct;

KKT_Residual_SOC.C_Fsb = KKT_Residual.C_Fsb + FunEval_full.C;
KKT_Residual_SOC.F_Fsb = KKT_Residual.F_Fsb + FunEval_full.F;

if (strcmp(OCPEC.VI_mode,'Reg_NCPs')) || (strcmp(OCPEC.VI_mode, 'Reg_Scholtes'))
    PHI_Fsb_Correct = solver.computeFB_minusInvPSIbPSI(FunEval.PSIphiPHI, FunEval_full.PSIphi);
    KKT_Residual_SOC.PHI_Fsb = KKT_Residual.PHI_Fsb + PHI_Fsb_Correct;
elseif strcmp(OCPEC.VI_mode,  'SmoothingEquation')
    KKT_Residual_SOC.PHI_Fsb = KKT_Residual.PHI_Fsb + FunEval_full.PHI;
end 

% compute second order correction direction and iterate
[dY_SOC, ~] = solver.SearchDirection_Riccati(KKT_Residual_SOC, KKT_Matrix);

Iterate_SOC.sigma  = Iterate.sigma  + dY_SOC(              1 : Dim.Node(1), :);
Iterate_SOC.eta    = Iterate.eta    + dY_SOC(Dim.Node(1) + 1 : Dim.Node(2), :);
Iterate_SOC.lambda = Iterate.lambda + dY_SOC(Dim.Node(2) + 1 : Dim.Node(3), :);
Iterate_SOC.gamma  = Iterate.gamma  + dY_SOC(Dim.Node(3) + 1 : Dim.Node(4), :);
Iterate_SOC.tau    = Iterate.tau    + dY_SOC(Dim.Node(4) + 1 : Dim.Node(5), :);
Iterate_SOC.x      = Iterate.x      + dY_SOC(Dim.Node(5) + 1 : Dim.Node(6), :);
Iterate_SOC.p      = Iterate.p      + dY_SOC(Dim.Node(6) + 1 : Dim.Node(7), :);
Iterate_SOC.w      = Iterate.w      + dY_SOC(Dim.Node(7) + 1 : Dim.Node(8), :);

end