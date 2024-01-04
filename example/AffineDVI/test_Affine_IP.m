clear all
clc

%% construct OCPEC problem
OCPEC = OCPEC_AffineDVI();

%% discretize OCPEC into a NLP problem
NLP = NLP_Formulation(OCPEC);

%% create solver
% problem
Prob = struct('x', NLP.z, 'f', NLP.J, 'g', [NLP.h; NLP.c; NLP.g], 'p', NLP.s);

% option (solver.print_options)
Option = struct;
Option.print_time = false;
Option.ipopt.max_iter = 2000;
Option.ipopt.tol = 1e-4;
Option.ipopt.print_level = 0;
Option.ipopt.mu_target = 0.5 * (1e-3)^2;

% create solver
solver = casadi.nlpsol('solver', 'ipopt', Prob, Option);

%% problem solve
z_Init = ones(NLP.Dim.z, 1);
s_Init = 1e-2;
s_End = 1e-8; 

[z_Opt, Info] = solveNLP_IP_homotopy(OCPEC, NLP, solver, z_Init, s_Init, s_End);

%%
plotResult_AffineDVI(OCPEC, NLP, z_Opt)