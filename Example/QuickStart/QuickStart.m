%%  A demo for quick start of NIPOCPEC using an affine DVI example
clear all
clc
delete Gen_InitialGuess.mat
[~, ~, ~] = rmdir('autoGen_CodeFiles', 's');

%% create an dynamical system
% specify dynamics variable dimension and using Class DifferentialVariationalInequalities to create an dynamical system
tau_Dim = 1;
x_Dim = 2;
p_Dim = 1;
plant = DifferentialVariationalInequalities(tau_Dim, x_Dim, p_Dim);

% set dynamics variable limit
tau_Max = 2;
tau_Min = -2;
x_Max = [2; 2];
x_Min = [-2; -2];
plant.setDynVarLimit(tau_Max, tau_Min, x_Max, x_Min);

% set state equation function
plant.f = [1, -3; -8, 10] * plant.x + [-3; -1] * plant.p + [4; 8] * plant.tau;

% set equilibrium constraint
plant.K = [1, -3] * plant.x + 5 * plant.p + 3 * plant.tau;
plant.l = -1;
plant.u = 1;

% set time step for discretization 
timeStep = 0.01;
plant.timeStep = timeStep;

% generate necessary code files
plant.codeGen();

%% formulate an OCPEC problem 
% set discretization
nStages = 100; 

% set initial state 
InitState = [-1/2; -1];

% set reference state and weight matrix in cost function
StageCost.xRef = repmat([0; 0], 1, nStages);
StageCost.tauRef = zeros(1, nStages);
StageCost.xWeight = [20; 20];
StageCost.tauWeight = 1; 
TerminalCost.xRef = [0; 0];
TerminalCost.tauRef = 0;
TerminalCost.xWeight = [20; 20];
TerminalCost.tauWeight = 1;

% specify the type of VI perturbed reformulation
VI_mode = 'Reg_Scholtes'; % option: 'SmoothingEquation', 'Reg_NCPs', 'Reg_Scholtes'

% using Class OCPEC_Perturbed to formulate an OCPEC problem 
OCPEC = OCPEC_Perturbed(plant, timeStep, nStages, InitState, StageCost, TerminalCost, VI_mode);

% set equality and inequality path constraint (optional)
C = [];
G = [];
OCPEC.setEqualityConstraints(C);
OCPEC.setInequalityConstraints(G);

% generate necessary code files
OCPEC.codeGen()

% show information
OCPEC.showInfo()

%% create a NIPOCPEC_Solver object 
% using Class NIPOCPEC_Solver to create a NIPOCPEC_Solver object
solver = NIPOCPEC_Solver(OCPEC);

% set solver option 
solver.Option.maxIterNum = 500;
solver.Option.Tolerance.KKT_Error_Total = 1e-2;
solver.Option.Tolerance.KKT_Error_Feasibility = 1e-4;
solver.Option.Tolerance.KKT_Error_Stationarity = 1e-4;
solver.Option.zInit = 1e-1; 
solver.Option.zEnd  = 1e-3;
solver.Option.sInit = 1e-1;
solver.Option.sEnd  = 1e-3;

% show solver information
solver.showInfo();

% generate necessary code files
solver.codeGen();

% generate initial guess
solver.generateInitialGuess();

%% solving OCPEC
% load initial guess
Gen_InitialGuess = load('Gen_InitialGuess.mat');
IterateInit = Gen_InitialGuess.Iterate;

% solve OCPEC
[solution, Info] = solver.solveOCPEC(IterateInit);

%% show result
% iteration information
solver.showResult(Info)

% plot trajectory
timeAxis = 0 : timeStep : nStages * timeStep;

K = zeros(plant.Dim.p, nStages);
for n = 1 : nStages
    K(:, n) = plant.computeVI_Function(solution.tau(:, n), solution.x(:, n), solution.p(:, n));
end

figure(111)
subplot(3,1,1)
plot(timeAxis, [InitState(1), solution.x(1, :)], 'r',...
     timeAxis, [InitState(2), solution.x(2, :)], 'g', 'LineWidth',1.2)
legend('x1', 'x2')
xlabel('time(s)')
title('system state')

subplot(3,1,2)
plot(timeAxis(2:end), solution.tau(1,:), 'LineWidth', 1.2)
xlabel('time(s)')
title('control')

subplot(3,1,3)
plot(timeAxis(2:end), solution.p(1, :), 'k',...
     timeAxis(2:end), K(1, :), 'b', 'LineWidth', 1.2)
legend('p', 'K') 
xlabel('time(s)')
title('equilibrium dynamics')
