clear all
clc
delete Acrobot.mp4

% parameter tracking is not satisfactory, problem in Hessian approximation due to the highly nonlinear dynamics?

%% create OCPEC, NLP and solver
OCPEC = OCPEC_Acrobot();
OCPEC.timeHorizon = 5;
OCPEC.nStages = 250;
OCPEC.timeStep = OCPEC.timeHorizon ./ OCPEC.nStages;
NLP = NLP_Formulation(OCPEC);

s_Init = 5e-2;
s_End = 1e-2;
sigma_Init = 1e-3;
sigma_End = 1e-4;

Option = NIP_C_Solver.createSolverOption();
Option.maxIterNum = 1000;
Option.KKT_scaling_max = 1;
Option.tol.KKT_error_primal = 1e-6; 
Option.tol.KKT_error_dual = 1e-2; 
Option.tol.KKT_error_complementarity = (sigma_Init)^2;
Option.tol.KKT_error_total = 1e-6; 
Option.tol.dYNorm = 1e-8;
Option.RegParam.nu_h = 1e-7; % rank-deficiency of equality constraint h
Option.RegParam.nu_c = 1e-7; % non-negative definite of diagonal matrix related to inequality constraint c
Option.RegParam.nu_g = 1e-7; % non-negative definite of diagonal matrix related to inequality constraint g
Option.RegParam.nu_H = 1e-6; % non-positive definite of Hessian matrix
Option.KKT.Hessian_approximation = 'Gauss_Newton'; % 'Exact', 'Gauss_Newton'
Option.Continuation.s_Init = s_Init;
Option.Continuation.s_End = s_End;
Option.Continuation.sigma_Init = sigma_Init;
Option.Continuation.sigma_End = sigma_End;
Option.Continuation.kappa_times = 0.95;
Option.Continuation.kappa_exp = 1.05;
Option.Continuation.tol.KKT_error = 1e-16;
Option.Continuation.tol.VI_nat_res = 1e-16;
Option.Continuation.AdditionNewtonStep = 1;

solver = NIP_C_Solver(OCPEC, NLP, Option);

%% solve problem
Z_Init = zeros(NLP.Dim.z_Node(4), OCPEC.nStages);
Z_Init(NLP.Dim.z_Node(1) + 1 : NLP.Dim.z_Node(2), :) = randn(OCPEC.Dim.u, OCPEC.nStages);
z_Init = reshape(Z_Init, [], 1);

[z_Opt, Info] = solver.solveNLP(z_Init);

%% show result
plotResult_Acrobot(OCPEC, NLP, z_Opt)

animateTrajectory_Acrobot(OCPEC, NLP, z_Opt)
