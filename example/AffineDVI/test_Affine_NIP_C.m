clear all
clc

%% create OCPEC, NLP and solver
OCPEC = OCPEC_AffineDVI();

NLP = NLP_Formulation(OCPEC);


s_Init = 1e-2;
s_End = 1e-8;
sigma_Init = 1e-2;
sigma_End = 1e-6;

Option = NIP_C_Solver.createSolverOption();
Option.maxIterNum = 1000;
Option.KKT_scaling_max = 1;
Option.tol.KKT_error_primal = 1e-6; 
Option.tol.KKT_error_dual = 1e-6; 
Option.tol.KKT_error_complementarity = (sigma_Init)^2;
Option.tol.KKT_error_total = 1e-4; 
Option.tol.dYNorm = 1e-8;
Option.RegParam.nu_h = 1e-6; % rank-deficiency of equality constraint h
Option.RegParam.nu_c = 1e-6; % non-negative definite of diagonal matrix related to inequality constraint c
Option.RegParam.nu_g = 1e-6; % non-negative definite of diagonal matrix related to inequality constraint g
Option.RegParam.nu_H = 1e-6; % non-positive definite of Hessian matrix
Option.KKT.Hessian_approximation = 'Gauss_Newton'; % 'Exact', 'Gauss_Newton'
Option.Continuation.s_Init = s_Init;
Option.Continuation.s_End = s_End;
Option.Continuation.sigma_Init = sigma_Init;
Option.Continuation.sigma_End = sigma_End;
Option.Continuation.kappa_times = 0.9;
Option.Continuation.kappa_exp = 1.1;
Option.Continuation.tol.KKT_error = 1e-16;
Option.Continuation.tol.VI_nat_res = 1e-16;
Option.Continuation.AdditionNewtonStep = 1;

solver = NIP_C_Solver(OCPEC, NLP, Option);

%% solve problem
z_Init = ones(NLP.Dim.z, 1);
[z_Opt, Info] = solver.solveNLP(z_Init);

%%
plotResult_AffineDVI(OCPEC, NLP, z_Opt)