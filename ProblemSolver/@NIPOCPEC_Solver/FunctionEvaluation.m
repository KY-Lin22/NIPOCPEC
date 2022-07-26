function FunEval = FunctionEvaluation(solver, Iterate, s, z, mode)
%FunctionEvaluation
%   Detailed explanation goes here

OCPEC = solver.OCPEC;

% init
FunEval = struct('L', [], 'G', [], 'C', [], 'F', [], 'PHI',[],...
    'PSIg', [], 'PSIgSigma', [], 'PSIgG', [], 'PSIphi', [], 'PSIphiGamma', [], 'PSIphiPHI', [],...
    'Lvar', [], 'Gvar', [], 'Cvar', [], 'Fvar', [], 'PHIvar',[],...
    'Hessian', []);
% cost function
FunEval.L = OCPEC.computeCostFun(Iterate, mode);

% constraint
[FunEval.G, FunEval.C, FunEval.F] = OCPEC.computeConstraintFun_G_C_F(Iterate);
FunEval.PHI = OCPEC.computeConstraintFun_PHI(Iterate, s);

% FB function and Jacobian for inequality-type constraint
[FunEval.PSIg, FunEval.PSIgSigma, FunEval.PSIgG] = solver.computeFB_Func_Jacobian(Iterate.sigma, FunEval.G, z);

if (strcmp(OCPEC.VI_mode,'Reg_NCPs')) || (strcmp(OCPEC.VI_mode, 'Reg_Scholtes'))
    [FunEval.PSIphi, FunEval.PSIphiGamma, FunEval.PSIphiPHI] = solver.computeFB_Func_Jacobian(Iterate.gamma, FunEval.PHI, z);
end

end

