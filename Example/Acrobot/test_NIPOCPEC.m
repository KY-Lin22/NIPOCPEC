%% Solving OCPEC with NIPOCPEC solver
clear all
clc
delete Acrobot.gif
delete Gen_InitialGuess.mat
[~, ~, ~] = rmdir('autoGen_CodeFiles', 's');
%% create an dynamics system 
timeStep = 0.05;
mass = [1; 1];
linkLength = [1; 1];
linkCenter = [0.5; 0.5];
inertia = [0.33; 0.33];
q2_min = -3/6*pi;
q2_max = 3/6*pi;
plant = Acrobot(timeStep, mass, linkLength, linkCenter, inertia, q2_min, q2_max);
tau_Max = 100;
tau_Min = -100;
x_Max = [2*pi; q2_max; 30; 30];
x_Min = [-2*pi; q2_min; -30; -30];
plant.setDynVarLimit(tau_Max, tau_Min, x_Max, x_Min) ;
plant.codeGen();

%% formulate an OCPEC problem 
nStages = 100;
InitState = [-0/3*pi; -0/3*pi; 0; 0];
StageCost.xRef = repmat([pi; 0; 0; 0], 1, nStages);
StageCost.tauRef = zeros(1, nStages);
StageCost.xWeight = [1; 1; 1; 1];
StageCost.tauWeight = 1; 
TerminalCost.xRef = [pi; 0; 0; 0];
TerminalCost.tauRef = 0;
TerminalCost.xWeight = [10; 10; 1; 1];
TerminalCost.tauWeight = 1;
% show initial and reference configuration of given plant
plotConfiguration(plant, InitState, TerminalCost.xRef)

OCPEC = OCPEC_Perturbed(plant, timeStep, nStages, InitState, StageCost, TerminalCost, 'Reg_Scholtes');% 'SmoothingEquation', 'Reg_NCPs', 'Reg_Scholtes'

OCPEC.setEqualityConstraints([]);
OCPEC.setInequalityConstraints([]);

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
solver.Option.Tolerance.KKT_Error_Feasibility = 1e-2;
solver.Option.Tolerance.KKT_Error_Stationarity = 1e-2;

solver.Option.RegularParam.nu_J = 1e-7;
solver.Option.RegularParam.nu_G = 1e-7;
solver.Option.RegularParam.nu_H = 0;

solver.Option.LineSearch.stepSize_Min = 0.001;
solver.Option.employFeasibilityRestorationPhase = true;

solver.Option.zInit = 1e-1; 
solver.Option.zEnd  = 1e-3;
solver.Option.sInit = 1e-1;
solver.Option.sEnd  = 1e-3;

% show solver information
solver.showInfo();
% generate initial guess
solver.generateInitialGuess();

%% solving OCPEC (time test)
% load initial guess
Gen_InitialGuess = load('Gen_InitialGuess.mat');
IterateInit = Gen_InitialGuess.Iterate;
% 
timeTest_Num = 1;
timeTest_TimeElapsed = 0;
timeTest_IterNum = 0;
solver.Option.printLevel = 2;
for i = 1 : timeTest_Num
    [solution, Info] = solver.solveOCPEC(IterateInit);
    timeTest_TimeElapsed = timeTest_TimeElapsed + Info.iterProcess.Time.total;
    timeTest_IterNum = timeTest_IterNum + Info.iterProcess.iterNum;
end
disp(['timeTest_iterations: ', num2str(timeTest_IterNum), '; ',...
      'timeTest_TimeElapsed: ', num2str(timeTest_TimeElapsed,'%10.3f'), ' s; ',...
      'timeTest_TimeAverage: ', num2str(1000 * timeTest_TimeElapsed /timeTest_IterNum, '%10.2f'), ' ms/Iter' ]);

%% show result 
plant.plotSimuResult(timeStep, InitState, solution.tau, solution.x, solution.p)
plant.animateTrajectory(timeStep, InitState, solution.tau, solution.x, solution.p)
solver.showResult(Info)
