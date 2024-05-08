clear all
clc

%% create OCPEC and NLP
OCPEC = OCPEC_CartPoleWithFriction();
OCPEC.timeHorizon = 3;
OCPEC.nStages = 300;
OCPEC.timeStep = OCPEC.timeHorizon ./ OCPEC.nStages;

NLP = NLP_Formulation(OCPEC);

%% create solver
% problem
Prob = struct('x', NLP.z, 'f', NLP.J, 'g', [NLP.h; NLP.c; NLP.g], 'p', NLP.s);

% option (solver.print_options)
Option = struct;
Option.print_time = false;
Option.ipopt.max_iter = 2000;
Option.ipopt.print_level = 0;
Option.ipopt.tol = 1e-6;
% create solver
solver = casadi.nlpsol('solver', 'ipopt', Prob, Option);

%% problem solve
Z_Init = zeros(NLP.Dim.z_Node(4), OCPEC.nStages);
Z_Init(NLP.Dim.z_Node(1) + 1 : NLP.Dim.z_Node(2), :) = randn(OCPEC.Dim.u, OCPEC.nStages);
z_Init = reshape(Z_Init, [], 1);
s_Init = 5e-1;
s_End = 1e-8; 

[z_Opt, Info] = solveNLP_IP_homotopy(OCPEC, NLP, solver, z_Init, s_Init, s_End);

% rec
rec.s_Init = s_Init;
rec.s_End = s_End;
rec.z_Init = z_Init;
rec.z_Opt = z_Opt;
rec.Info = Info;
save('Data_IP.mat', 'rec')