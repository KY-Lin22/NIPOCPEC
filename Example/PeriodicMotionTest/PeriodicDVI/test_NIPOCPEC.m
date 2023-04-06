%% Solving OCPEC with NIPOCPEC solver
clear all
clc
delete holf_oscillator.gif
delete Gen_InitialGuess.mat
[~, ~, ~] = rmdir('autoGen_CodeFiles', 's');

%% create a dynamics system 
timeStep = 0.01;
OscillatorParam.convergeSpeed1 = 50; % control the speed for the oscillator to converge to the limit cycle (rho_1)
OscillatorParam.convergeSpeed2 = 50; % control the speed for the oscillator to converge to the limit cycle (rho_2)
OscillatorParam.radiusLimitCycle = 1; % radius of the limit cycle
OscillatorParam.periodLimitCycle = 0.3; % period of the limit cycle
OscillatorParam.dutyFactor = 0.4; % duty factor determining the fraction of rho_2 > 0 in a whole period
OscillatorParam.radiusSmoothVaries = 50; % ensure radiusLimitCycle varies smoothly across different half planes
OscillatorParam.InitPhase = [0; 0; pi; pi];
FootVelFuncParam.k_N = 0.1;
FootVelFuncParam.k_T = 0.1;
plant = PeriodicDVI(timeStep, OscillatorParam, FootVelFuncParam);

plant.codeGen();

%% formulate an OCPEC problem 
nStages = 60;
x_Init = 0.3;
linkLength1 = 0.267;
num_Osci = length(OscillatorParam.InitPhase);
InitState = zeros(4 * num_Osci, 1);
for i = 1 : num_Osci
    InitPhase_i = OscillatorParam.InitPhase(i);
    InitRadius = OscillatorParam.radiusLimitCycle;
    InitState_HolfOsc_i = [InitRadius * cos(InitPhase_i);...
        InitRadius * sin(InitPhase_i)];% holf oscillator
    % for the case that num_Osci = 4
    if (i == 1) || (i == 2)
        % front leg
        InitState_fpt_i = [0; x_Init + linkLength1];
    elseif (i == 3) || (i == 4)
        % hind leg
        InitState_fpt_i = [0; x_Init];
    end    
    InitState(1 + (i - 1) * 4 : 4 + (i - 1) * 4) = [InitState_HolfOsc_i; InitState_fpt_i];
end
[state, phase, radius] = plant.systemSimulation(InitState, nStages, timeStep);
TerminalState = state(:, end);
%%
plant.animateSystemSimulationResult(InitState, state, nStages, timeStep)

%%
timeAxis = 0 : timeStep : nStages * timeStep;
figure(1)
subplot(4,1,1)
plot(timeAxis, cos([OscillatorParam.InitPhase(1), phase(1,:)]), 'r', 'LineWidth',1.2)
legend('osci 1')
subplot(4,1,2)
plot(timeAxis, cos([OscillatorParam.InitPhase(2), phase(2,:)]), 'g', 'LineWidth',1.2)
legend('osci 2')
subplot(4,1,3)
plot(timeAxis, cos([OscillatorParam.InitPhase(3), phase(3,:)]), 'b', 'LineWidth',1.2)
legend('osci 3')
subplot(4,1,4)
plot(timeAxis, cos([OscillatorParam.InitPhase(4), phase(4,:)]), 'k', 'LineWidth',1.2)
legend('osci 4')
xlabel('time [s]')
%%
% cost function
StageCost.xRef = zeros(plant.Dim.x, nStages);
StageCost.tauRef = repmat(zeros(plant.Dim.tau, 1), 1, nStages);
StageCost.xWeight = 1 * ones(plant.Dim.x, 1);
StageCost.tauWeight = 10 * ones(plant.Dim.tau, 1); 

TerminalCost.xRef = TerminalState;
TerminalCost.tauRef = zeros(plant.Dim.tau, 1);
TerminalCost.xWeight = 1 * ones(plant.Dim.x, 1);
TerminalCost.tauWeight = 10 * ones(plant.Dim.tau, 1);

OCPEC = OCPEC_Perturbed(plant, timeStep, nStages, InitState, StageCost, TerminalCost, 'Reg_Scholtes');% 'SmoothingEquation', 'Reg_NCPs', 'Reg_Scholtes' 

C = [];
OCPEC.setEqualityConstraints(C);
G_gap = sym('G_gap', [num_Osci, 1]);
G_rho_gap_comple = sym('G_rho_gap_comple', [num_Osci, 1]);
for i = 1 : num_Osci
    G_gap(i) = plant.x(3 + (i - 1) * num_Osci);
    G_rho_gap_comple(i) = 0.1 - plant.x(2 + (i - 1) * num_Osci) * plant.x(3 + (i - 1) * num_Osci);
end
G = [G_gap;...
    G_rho_gap_comple];
OCPEC.setInequalityConstraints(G);

% generate necessary code files
OCPEC.codeGen()

OCPEC.showInfo()
%% create a NIPOCPEC_Solver object 
solver = NIPOCPEC_Solver(OCPEC);
% generate necessary code files
solver.codeGen();

%% set option and generate initial guess
solver.Option.maxIterNum = 500;
solver.Option.Tolerance.KKT_Error_Total = 1e-2;
solver.Option.Tolerance.KKT_Error_Feasibility = 1e-4;
solver.Option.Tolerance.KKT_Error_Stationarity = 1e-2;

solver.Option.LineSearch.stepSize_Min = 0.001;
solver.Option.employFeasibilityRestorationPhase = true;
solver.Option.FRP.maxIterNum = 50;
solver.Option.zInit = 1e-1; 
solver.Option.zEnd  = 1e-2;
solver.Option.sInit = 1e-1;
solver.Option.sEnd  = 1e-2;

% show solver information
solver.showInfo();
% generate initial guess
solver.generateInitialGuess();

%% solving OCPEC (time test)
% load initial guess
Gen_InitialGuess = load('Gen_InitialGuess.mat');
IterateInit = Gen_InitialGuess.Iterate;
%IterateInit.tau = zeros(plant.Dim.tau, nStages);
IterateInit.x = state;
[solution, Info] = solver.solveOCPEC(IterateInit);

%%
plant.animateSystemSimulationResult(InitState, solution.x, nStages, timeStep)

figure(1)
for i = 1 : num_Osci
    subplot(num_Osci,1,i)
    footTraj_X_i = [InitState(4 + (i - 1)*4, :), solution.x(4 + (i - 1)*4, :)];
    footTraj_Y_i = [InitState(3 + (i - 1)*4, :), solution.x(3 + (i - 1)*4, :)];
    plot(footTraj_X_i,footTraj_Y_i, 'r', 'LineWidth',1.2)
    axisLimit_X = [x_Init - 0.1; x_Init + 4 * linkLength1 + 0.1];
    axisLimit_Y = [min(footTraj_Y_i); max(footTraj_Y_i) + 0.1];
    axis([axisLimit_X; axisLimit_Y]);
end

save('PeriodicDVI_InitGuess.mat', 'InitState', 'solution', 'OscillatorParam', 'FootVelFuncParam', 'x_Init');