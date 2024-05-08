clear all
clc

%% create solver set based on nStage sequence
timeHorizon = 3;
nStage_sequence = {60, 100, 300, 600};
name = {...
    '$\Delta t = 5 \cdot 10^{-2} (IP)$',...
    '$\Delta t = 3 \cdot 10^{-2} (IP)$',...
    '$\Delta t = 1 \cdot 10^{-2} (IP)$',...
    '$\Delta t = 5 \cdot 10^{-3} (IP)$'};
s_Init = 5e-1; 
s_End = 1e-8; 
OCPEC_set = cell(1, numel(nStage_sequence));
NLP_set = cell(1, numel(nStage_sequence));
solver_set = cell(1, numel(nStage_sequence));
for i = 1 : numel(nStage_sequence)
    % create OCPEC, NLP and solver
    OCPEC_i = OCPEC_CartPoleWithFriction();
    OCPEC_i.timeHorizon = timeHorizon;
    OCPEC_i.nStages = nStage_sequence{i};
    OCPEC_i.timeStep = OCPEC_i.timeHorizon ./ OCPEC_i.nStages;
    NLP_i = NLP_Formulation(OCPEC_i);
    Prob_i = struct('x', NLP_i.z, 'f', NLP_i.J, 'g', [NLP_i.h; NLP_i.c; NLP_i.g], 'p', NLP_i.s);
    Option_i = struct;
    Option_i.print_time = false;
    Option_i.ipopt.max_iter = 2000;
    Option_i.ipopt.print_level = 0;
    solver_i = casadi.nlpsol('solver', 'ipopt', Prob_i, Option_i);
    % save
    OCPEC_set{i} = OCPEC_i;
    NLP_set{i} = NLP_i;
    solver_set{i} = solver_i;
end

%% run test
% init record
rec.timeHorizon = timeHorizon;
rec.nStage_sequence = nStage_sequence;
rec.name = name;
rec.s_Init = s_Init;
rec.s_End = s_End;
rec.z_Init = cell(1, numel(nStage_sequence));
rec.z_Opt = cell(1, numel(nStage_sequence));
rec.Info = cell(1, numel(nStage_sequence));

% run
for i = 1 : numel(nStage_sequence)
    OCPEC_i = OCPEC_set{i};
    NLP_i = NLP_set{i};
    solver_i = solver_set{i};
    Z_Init_i = ones(NLP_i.Dim.z_Node(4), OCPEC_i.nStages);
    Z_Init_i(NLP_i.Dim.z_Node(1) + 1 : NLP_i.Dim.z_Node(2), :) = randn(OCPEC_i.Dim.u, OCPEC_i.nStages);
    z_Init_i = reshape(Z_Init_i, [], 1);

    [z_Opt_i, Info_i] = solveNLP_IP_homotopy(OCPEC_i, NLP_i, solver_i, z_Init_i, s_Init, s_End);
    rec.z_Init{i} = z_Init_i;
    rec.z_Opt{i} = z_Opt_i;
    rec.Info{i} = Info_i;
end
save('Data_time_step_IP.mat', 'rec')