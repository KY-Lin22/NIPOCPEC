clear all
clc
delete ThreeCarts.mp4

%% create OCPEC, NLP and solver
OCPEC = OCPEC_ThreeCartsSystem();
OCPEC.timeHorizon = 4;
OCPEC.nStages = 200;
OCPEC.timeStep = OCPEC.timeHorizon ./ OCPEC.nStages;
NLP = NLP_Formulation(OCPEC);

solver = NIP_C_Solver(OCPEC, NLP);

%% set option and solve problem
p_Init = [5e-1; 1e-1]; % [s; sigma]
p_End = [1e-8; 1e-6]; % [s; sigma]
solver.Option.maxIterNum = 1000;
solver.Option.KKT_scaling_max = 1;
solver.Option.tol.KKT_error_primal = 1e-6; 
solver.Option.tol.KKT_error_dual = 1e-2; 
solver.Option.tol.KKT_error_complementarity = (p_Init(2))^2;
solver.Option.tol.KKT_error_total = 1e-6; 
solver.Option.tol.dYNorm = 1e-8;
solver.Option.RegParam.nu_h = 1e-7; % rank-deficiency of equality constraint h
solver.Option.RegParam.nu_c = 1e-7; % non-negative definite of diagonal matrix related to inequality constraint c
solver.Option.RegParam.nu_g = 1e-7; % non-negative definite of diagonal matrix related to inequality constraint g
solver.Option.RegParam.nu_H = 1e-6; % non-positive definite of Hessian matrix
solver.Option.printLevel = 2;
solver.Option.NewtonCorrection.StepNum = 1;

% Euler-Newton
Z_Init = ones(NLP.Dim.z_Node(4), OCPEC.nStages);
Z_Init(NLP.Dim.z_Node(1) + 1 : NLP.Dim.z_Node(2), :) = randn(OCPEC.Dim.u, OCPEC.nStages);
z_Init = reshape(Z_Init, [], 1);

[z_Opt, Info] = solver.solveNLP(z_Init, p_Init, p_End);

%% show result
plotResult_ThreeCartsSystem(OCPEC, NLP, z_Opt)

animateTrajectory_ThreeCartsSystem(OCPEC, NLP, z_Opt)
