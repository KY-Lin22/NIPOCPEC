function KKT_Residual = computeKKT_Residual(solver, Iterate, FunEval)
%computeKKT_Residual
%   Detailed explanation goes here

OCPEC = solver.OCPEC;
Dim = OCPEC.Dim;
nStages = OCPEC.nStages;

% init
KKT_Residual = struct('G_Fsb', [], 'C_Fsb', [], 'F_Fsb', [], 'PHI_Fsb', [],...
    'HtauT', [], 'HxTlambdaNext', [], 'HpT', [], 'HwT', []);

%% primal feasibility
KKT_Residual.G_Fsb = solver.computeFB_minusInvPSIbPSI(FunEval.PSIgG, FunEval.PSIg);

KKT_Residual.C_Fsb = FunEval.C;

KKT_Residual.F_Fsb = FunEval.F;

if (strcmp(OCPEC.VI_mode,'Reg_NCPs')) || (strcmp(OCPEC.VI_mode, 'Reg_Scholtes'))
    KKT_Residual.PHI_Fsb = solver.computeFB_minusInvPSIbPSI(FunEval.PSIphiPHI, FunEval.PSIphi);
elseif strcmp(OCPEC.VI_mode,  'SmoothingEquation')
    KKT_Residual.PHI_Fsb = FunEval.PHI;
end 

%% dual feasibility
lambdaNext = [Iterate.lambda(:, 2 : end), zeros(Dim.lambda, 1)];

KKT_Residual.HtauT         = zeros(Dim.tau, nStages);
KKT_Residual.HxTlambdaNext = zeros(Dim.x, nStages);
KKT_Residual.HpT           = zeros(Dim.p, nStages);
KKT_Residual.HwT           = zeros(Dim.w, nStages);
for n = 1 : nStages
    % load
    sigma_n = Iterate.sigma(:, n);
    eta_n = Iterate.eta(:, n);
    lambda_n = Iterate.lambda(:, n);
    gamma_n = Iterate.gamma(:, n);
    
    Ltau_n = FunEval.Lvar.Ltau(:, (n - 1) * Dim.tau + 1: n * Dim.tau);
    Lx_n = FunEval.Lvar.Lx(:, (n - 1) * Dim.x + 1 : n * Dim.x);
    Lp_n = FunEval.Lvar.Lp(:, (n - 1) * Dim.p + 1 : n * Dim.p);
    Lw_n = FunEval.Lvar.Lw(:, (n - 1) * Dim.w + 1 : n * Dim.w);
    
    Gtau_n = FunEval.Gvar.Gtau(:, (n - 1) * Dim.tau + 1 : n * Dim.tau);
    Gx_n = FunEval.Gvar.Gx(:, (n - 1) * Dim.x + 1 : n * Dim.x);
    Gp_n = FunEval.Gvar.Gp(:, (n - 1) * Dim.p + 1 : n * Dim.p);
    
    Ctau_n = FunEval.Cvar.Ctau(:, (n - 1) * Dim.tau + 1 : n * Dim.tau);
    Cx_n = FunEval.Cvar.Cx(:, (n - 1) * Dim.x + 1 : n * Dim.x);
    Cp_n = FunEval.Cvar.Cp(:, (n - 1) * Dim.p + 1 : n * Dim.p);
    Cw_n = FunEval.Cvar.Cw(:, (n - 1) * Dim.w + 1 : n * Dim.w);
    
    Ftau_n = FunEval.Fvar.Ftau(:, (n - 1) * Dim.tau + 1 : n * Dim.tau);
    Fx_n = FunEval.Fvar.Fx(:, (n - 1) * Dim.x + 1 : n * Dim.x);
    Fp_n = FunEval.Fvar.Fp(:, (n - 1) * Dim.p + 1 : n * Dim.p);
    
    PHIp_n = FunEval.PHIvar.PHIp(:, (n - 1) * Dim.p + 1 : n * Dim.p);
    PHIw_n = FunEval.PHIvar.PHIw(:, (n - 1) * Dim.w + 1 : n * Dim.w);   
    
    % compute dual feasibility
    if (strcmp(OCPEC.VI_mode, 'Reg_Scholtes')) || (strcmp(OCPEC.VI_mode, 'Reg_NCPs'))
        Htau_n = Ltau_n - sigma_n' * Gtau_n + eta_n' * Ctau_n + lambda_n' * Ftau_n;
        Hx_n   = Lx_n   - sigma_n' * Gx_n   + eta_n' * Cx_n   + lambda_n' * Fx_n;
        Hp_n   = Lp_n   - sigma_n' * Gp_n   + eta_n' * Cp_n   + lambda_n' * Fp_n  - gamma_n' * PHIp_n;
        Hw_n   = Lw_n                       + eta_n' * Cw_n                       - gamma_n' * PHIw_n;
    elseif strcmp(OCPEC.VI_mode, 'SmoothingEquation')
        Htau_n = Ltau_n - sigma_n' * Gtau_n + eta_n' * Ctau_n + lambda_n' * Ftau_n;
        Hx_n   = Lx_n   - sigma_n' * Gx_n   + eta_n' * Cx_n   + lambda_n' * Fx_n;
        Hp_n   = Lp_n   - sigma_n' * Gp_n   + eta_n' * Cp_n   + lambda_n' * Fp_n  + gamma_n' * PHIp_n;
        Hw_n   = Lw_n                       + eta_n' * Cw_n                       + gamma_n' * PHIw_n;
    end     
    
    % record  
    KKT_Residual.HtauT(:, n) = Htau_n';
    KKT_Residual.HxTlambdaNext(:, n) = Hx_n' + lambdaNext(:, n);
    KKT_Residual.HpT(:, n) = Hp_n';   
    KKT_Residual.HwT(:, n) = Hw_n';        
end

end

