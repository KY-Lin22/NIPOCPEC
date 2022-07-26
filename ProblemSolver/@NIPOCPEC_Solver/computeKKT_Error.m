function KKT_Error = computeKKT_Error(solver, Iterate, FunEval, KKT_Residual)
%computeKKT_Error
%   Detailed explanation goes here

OCPEC = solver.OCPEC;

% compute scaling parameter for KKT stationarity
scaling_max = 100;
LAMBDA_norm = norm(reshape(Iterate.sigma, [], 1), 1) + norm(reshape(Iterate.eta, [], 1), 1) ...
             + norm(reshape(Iterate.lambda, [], 1), 1) + norm(reshape(Iterate.gamma, [], 1), 1);          
LAMBDA_trial = (LAMBDA_norm)/(OCPEC.nStages * OCPEC.Dim.LAMBDA);
stationarityScaling = max([scaling_max, LAMBDA_trial])/scaling_max;

% compute KKT_Error
if (strcmp(OCPEC.VI_mode,'Reg_NCPs')) || (strcmp(OCPEC.VI_mode, 'Reg_Scholtes'))
    Feasibility = [FunEval.PSIg; FunEval.C; FunEval.F; FunEval.PSIphi];
elseif strcmp(OCPEC.VI_mode,  'SmoothingEquation')
    Feasibility = [FunEval.PSIg; FunEval.C; FunEval.F; FunEval.PHI];
end 

Stationarity = [KKT_Residual.HtauT; KKT_Residual.HxTlambdaNext; KKT_Residual.HpT; KKT_Residual.HwT];

KKT_Error.Feasibility = norm(reshape(Feasibility, [], 1), Inf);

KKT_Error.Stationarity = norm(reshape(Stationarity, [], 1), Inf)/stationarityScaling;

KKT_Error.Total = max([KKT_Error.Feasibility, KKT_Error.Stationarity]);
end

