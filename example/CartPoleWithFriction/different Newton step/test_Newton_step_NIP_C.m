clear all
clc
%% create solver set based on nStage sequence
timeHorizon = 3;
nStages = 300;
num_Newton_step = {1, 2, 3};
name = {'NIP (1 corrector step)', 'NIP (2 corrector step)', 'NIP (3 corrector step)'};
p_Init = [5e-1; 1e-1]; % [s; sigma]
p_End = [1e-8; 1e-6]; % [s; sigma]
solver_set = cell(1, numel(num_Newton_step));
for i = 1 : numel(num_Newton_step)
    % create OCPEC, NLP and solver
    OCPEC_i = OCPEC_CartPoleWithFriction();
    OCPEC_i.timeHorizon = timeHorizon;
    OCPEC_i.nStages = nStages;
    OCPEC_i.timeStep = OCPEC_i.timeHorizon ./ OCPEC_i.nStages;
    NLP_i = NLP_Formulation(OCPEC_i);
    solver_i = NIP_C_Solver(OCPEC_i, NLP_i);
    % solver option
    solver_i.Option.maxIterNum = 1000;
    solver_i.Option.KKT_scaling_max = 1;
    solver_i.Option.tol.KKT_error_primal = 1e-6;
    solver_i.Option.tol.KKT_error_dual = 1e-2;
    solver_i.Option.tol.KKT_error_complementarity = (p_Init(2))^2;
    solver_i.Option.tol.KKT_error_total = 1e-6;
    solver_i.Option.tol.dYNorm = 1e-8;
    solver_i.Option.RegParam.nu_h = 1e-7; % rank-deficiency of equality constraint h
    solver_i.Option.RegParam.nu_c = 1e-7; % non-negative definite of diagonal matrix related to inequality constraint c
    solver_i.Option.RegParam.nu_g = 1e-7; % non-negative definite of diagonal matrix related to inequality constraint g
    solver_i.Option.RegParam.nu_H = 1e-6; % non-positive definite of Hessian matrix
    solver_i.Option.printLevel = 2;
    solver_i.Option.NewtonCorrection.StepNum = num_Newton_step{i};
    % save
    solver_set{i} = solver_i;
end
%% run test
% init record
rec.timeHorizon = timeHorizon;
rec.nStages = nStages;
rec.num_Newton_step = num_Newton_step;
rec.name = name;
rec.p_Init = p_Init;
rec.p_End = p_End;
rec.z_Init = cell(1, numel(num_Newton_step));
rec.z_Opt = cell(1, numel(num_Newton_step));
rec.Info = cell(1, numel(num_Newton_step));
% run
for i = 1 : numel(num_Newton_step)
    solver_i = solver_set{i};
    Z_Init_i = ones(solver_i.NLP.Dim.z_Node(4), solver_i.OCPEC.nStages);
    Z_Init_i(solver_i.NLP.Dim.z_Node(1) + 1 : solver_i.NLP.Dim.z_Node(2), :) = randn(solver_i.OCPEC.Dim.u, solver_i.OCPEC.nStages);
    z_Init_i = reshape(Z_Init_i, [], 1);

    [z_Opt_i, Info_i] = solver_i.solveNLP(z_Init_i, p_Init, p_End);
    rec.z_Init{i} = z_Init_i;
    rec.z_Opt{i} = z_Opt_i;
    rec.Info{i} = Info_i;
end
save('Data_Newton_step_NIP.mat', 'rec')