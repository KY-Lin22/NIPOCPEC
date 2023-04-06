%% Solving OCPEC with NIPOCPEC solver
clear all
clc
delete Gen_InitialGuess.mat
[~, ~, ~] = rmdir('autoGen_CodeFiles', 's');

%% create an dynamics system
timeStep = 0.02;
plant = FilippovDI_BVI(timeStep);
plant.codeGen();

%% formulate an OCPEC problem 
nStages = 100;
InitState = -1;

StageCost.xRef = zeros(1, nStages);
StageCost.tauRef = zeros(plant.Dim.tau, nStages);
StageCost.xWeight = 2;
StageCost.tauWeight = 0.00001*zeros(plant.Dim.tau, 1); 
TerminalCost.xRef = 5/3;
TerminalCost.tauRef = zeros(plant.Dim.tau, 1);
TerminalCost.xWeight = 2;
TerminalCost.tauWeight = 0.00001*zeros(plant.Dim.tau, 1);

OCPEC = OCPEC_Perturbed(plant, timeStep, nStages, InitState, StageCost, TerminalCost, 'Reg_Scholtes');% 'SmoothingEquation', 'Reg_NCPs', 'Reg_Scholtes' 

% set C
OCPEC.setEqualityConstraints([]);

% set G
OCPEC.setInequalityConstraints([]);

% generate necessary code files
OCPEC.codeGen()

OCPEC.showInfo()

%% create a NIPOCPEC_Solver object 
solver = NIPOCPEC_Solver(OCPEC);
% generate necessary code files
solver.codeGen();

%% set option and generate initial guess
solver.Option.maxIterNum = 100;
solver.Option.Tolerance.KKT_Error_Total = 1e-2;
solver.Option.Tolerance.KKT_Error_Feasibility = 1e-4;
solver.Option.Tolerance.KKT_Error_Stationarity = 1e-4;

solver.Option.LineSearch.stepSize_Min = 0.01;
solver.Option.employFeasibilityRestorationPhase = true;

solver.Option.zInit = 1e-1; 
solver.Option.zEnd  = 1e-4;
solver.Option.sInit = 1e-1;
solver.Option.sEnd  = 1e-4;

% show solver information
solver.showInfo();
% generate initial guess
solver.generateInitialGuess();

%% solving OCPEC (time test)
% load initial guess
Gen_InitialGuess = load('Gen_InitialGuess.mat');
IterateInit = Gen_InitialGuess.Iterate;

[solution, Info] = solver.solveOCPEC(IterateInit);

%% show result 
plant.plotSimuResult(timeStep, InitState, solution.tau, solution.x, solution.p)
