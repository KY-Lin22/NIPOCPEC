clear all
clc

%% construct OCPEC problem
OCPEC = OCPEC_CartPoleWithFriction();

%% discretize OCPEC into a NLP problem
NLP = NLP_Formulation(OCPEC);

%% create solver
solver = NIP_C_Solver(OCPEC, NLP);

%% set option and solve problem
% one solve
solver.Option.maxIterNum = 1000;
solver.Option.KKT_scaling_max = 1;
solver.Option.tol.KKT_error_primal = 1e-6; 
solver.Option.tol.KKT_error_dual = 1e-2; 
solver.Option.tol.KKT_error_total = 1e-6; 
solver.Option.tol.dYNorm = 1e-8;
solver.Option.RegParam.nu_h = 1e-7; % rank-deficiency of equality constraint h
solver.Option.RegParam.nu_c = 1e-7; % non-negative definite of diagonal matrix related to inequality constraint c
solver.Option.RegParam.nu_g = 1e-7; % non-negative definite of diagonal matrix related to inequality constraint g
solver.Option.RegParam.nu_H = 1e-6; % non-positive definite of Hessian matrix
solver.Option.recordLevel = 1;
solver.Option.printLevel = 2;
solver.Option.NewtonCorrection.StepNum = 1;
% single step
% z_Init = ones(NLP.Dim.z, 1);
% p = [1e-8; 1e-3]; % [s; sigma]
% [z_Opt, Info] = solver.non_interior_point_method(z_Init, p);

% Euler-Newton
Z_Init = ones(NLP.Dim.z_Node(4), OCPEC.nStages);
Z_Init(NLP.Dim.z_Node(1) + 1 : NLP.Dim.z_Node(2), :) = randn(OCPEC.Dim.u, OCPEC.nStages);
z_Init = reshape(Z_Init, [], 1);
% save('initial_guess.mat', 'z_Init');
% initial_guess = load('initial_guess.mat');
% z_Init = initial_guess.z_Init;
p_Init = [5e-1; 1e-1]; % [s; sigma]
p_End = [1e-8; 1e-6]; % [s; sigma]
[z_Opt, Info] = solver.solveNLP(z_Init, p_Init, p_End);
%%
% save('Data_T_3s_t_100ms_NIP.mat', 'z_Opt', 'Info', 'z_Init', 'p_Init', 'p_End')
save('Data_corrector_3step_NIP.mat', 'z_Opt', 'Info', 'z_Init', 'p_Init', 'p_End')

%%
solver.showResult(Info)

plotResult_CartPoleWithFriction(OCPEC, NLP, z_Opt)
