clear all
clc

%% construct OCPEC problem
OCPEC = OCPEC_AffineDVI();

%% discretize OCPEC into a NLP problem
NLP = NLP_Formulation(OCPEC);

%% create solver
solver = NIP_C_Solver(OCPEC, NLP);

%% set option and solve problem
% one solve
solver.Option.KKT_scaling_max = 100;
solver.Option.tol.KKT_error_primal = 1e-6; 
solver.Option.tol.KKT_error_dual = 1e-6; 
solver.Option.tol.KKT_error_total = 1e-4; 
solver.Option.tol.dYNorm = 1e-8;
solver.Option.RegParam.nu_h = 1e-6; % rank-deficiency of equality constraint h
solver.Option.RegParam.nu_c = 1e-6; % non-negative definite of diagonal matrix related to inequality constraint c
solver.Option.RegParam.nu_g = 1e-6; % non-negative definite of diagonal matrix related to inequality constraint g
solver.Option.RegParam.nu_H = 1e-6; % non-positive definite of Hessian matrix
solver.Option.recordLevel = 1;
solver.Option.printLevel = 2;
solver.Option.NewtonCorrection.StepNum = 2;
% single step
% z_Init = ones(NLP.Dim.z, 1);
% p = [1e-8; 1e-3]; % [s; sigma]
% [z_Opt, Info] = solver.non_interior_point_method(z_Init, p);

% Euler-Newton
z_Init = ones(NLP.Dim.z, 1);
p_Init = [1e-2; 1e-2]; % [s; sigma]
p_End = [1e-8; 1e-3]; % [s; sigma]
[z_Opt, Info] = solver.solveNLP(z_Init, p_Init, p_End);

%%
solver.showResult(Info)

plotResult_AffineDVI(OCPEC, NLP, z_Opt)