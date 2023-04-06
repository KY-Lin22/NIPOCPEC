%% Solving OCPEC with NIPOCPEC solver
clear all
clc
delete Hopper.gif
delete Gen_InitialGuess.mat
[~, ~, ~] = rmdir('autoGen_CodeFiles', 's');

%% create an dynamics system 
timeStep = 0.01;
frictionCoeff = 0.45;
BodyParam.mass = [1; 0.1];
BodyParam.inertia = [0.25; 0.025];
OscillatorParam.convergeSpeed1 = 50; % control the speed for the oscillator to converge to the limit cycle (rho_1)
OscillatorParam.convergeSpeed2 = 50; % control the speed for the oscillator to converge to the limit cycle (rho_2)
OscillatorParam.radiusLimitCycle = 1; % radius of the limit cycle
OscillatorParam.periodLimitCycle = 0.6; % period of the limit cycle
OscillatorParam.dutyFactor = 0.5; % duty factor determining the fraction of rho_2 > 0 in a whole period
OscillatorParam.radiusSmoothVaries = 50; % ensure radiusLimitCycle varies smoothly across different half planes

plant = Hopper_PeriodicGait(timeStep, frictionCoeff, BodyParam, OscillatorParam);

tau_Max = [50; 50; 100];
tau_Min = [-50; -50; -100];
x_Max = [0.8; 1.5; pi; 0.5; 10; 10; 5; 5; 100; 100; 100; 100];
x_Min = [0; 0; -pi; 0.1; -10; -10; -5; -5; -100; -100; -100; -100];
plant.setDynVarLimit(tau_Max, tau_Min, x_Max, x_Min) ;
plant.codeGen();

%% formulate an OCPEC problem 
nStages = 120;

InitState_rbg = [0.1; 0.5; 0; 0.5; 0; 0; 0; 0];
InitPhase = 0;
InitRadius = OscillatorParam.radiusLimitCycle;
InitState_HolfOsc = [InitRadius * cos(InitPhase);...
    InitRadius * sin(InitPhase)];% holf oscillator
InitState_fpt = [InitState_rbg(2) - InitState_rbg(4) * cos(InitState_rbg(3));...
    InitState_rbg(1) + InitState_rbg(4) * sin(InitState_rbg(3))]; % foot position

InitState = [InitState_rbg;...
    InitState_HolfOsc;...
    InitState_fpt];

RefState_rbg = [0.7; 0.5; 0; 0.5; 0; 0; 0; 0];
RefPhase = 0;
RefRadius = OscillatorParam.radiusLimitCycle;
RefState_HolfOsc = [RefRadius * cos(RefPhase);...
    RefRadius * sin(RefPhase)];% holf oscillator
RefState_fpt = [RefState_rbg(2) - RefState_rbg(4) * cos(RefState_rbg(3));...
    RefState_rbg(1) + RefState_rbg(4) * sin(RefState_rbg(3))]; % foot position
RefState = [RefState_rbg;...
    RefState_HolfOsc;...
    RefState_fpt];

xRef = TrajectoryInterpolation(InitState, RefState, nStages);
StageCost.xRef = xRef;
StageCost.tauRef = repmat([0; 0; 0], 1, nStages);
StageCost.xWeight = [100; 1; 1; 1; 0.1; 0.1; 0.1; 0.1;...
    0.1; 0.1; 1; 10];
StageCost.tauWeight = [0.1; 0.1; 0.001]; 

TerminalCost.xRef = RefState;
TerminalCost.tauRef = [0; 0; 0];
TerminalCost.xWeight = [100; 1; 1; 1; 0.1; 0.1; 0.1; 0.1;...
    0.1; 0.1; 1; 10];
TerminalCost.tauWeight = [0.1; 0.1; 0.001];

OCPEC = OCPEC_Perturbed(plant, timeStep, nStages, InitState, StageCost, TerminalCost, 'Reg_Scholtes');% 'SmoothingEquation', 'Reg_NCPs', 'Reg_Scholtes'  

W_T = [1, 0, OCPEC.x(4)*cos(OCPEC.x(3)), sin(OCPEC.x(3))];
vel_T = W_T * OCPEC.x(5:8);
C1 = OCPEC.p(2) - OCPEC.p(3) - vel_T; % using two auxilary variable for vel_T to reformulate friction 

C2 = OCPEC.x(11) - (OCPEC.x(2) - OCPEC.x(4) * cos(OCPEC.x(3)));% relationship between foot position variable and rigid body variable
C3 = OCPEC.x(12) - (OCPEC.x(1) + OCPEC.x(4) * sin(OCPEC.x(3)));

C = [C1; C2; C3];
OCPEC.setEqualityConstraints(C);

G = 0.01 - OCPEC.p(1) * vel_T; % penalty slip motion
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
solver.Option.Tolerance.KKT_Error_Total = 1e-1;
solver.Option.Tolerance.KKT_Error_Feasibility = 1e-2;
solver.Option.Tolerance.KKT_Error_Stationarity = 1e-2;

solver.Option.RegularParam.nu_J = 1e-7;
solver.Option.RegularParam.nu_G = 1e-7;
solver.Option.RegularParam.nu_H = 0;

solver.Option.LineSearch.stepSize_Min = 0.005;
solver.Option.employFeasibilityRestorationPhase = true;
solver.Option.FRP.maxIterNum = 50;
solver.Option.FRP.betaInit = 1;
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

[solution, Info] = solver.solveOCPEC(IterateInit);

%% show result 
plant.plotSimuResult(timeStep, InitState, solution.tau, solution.x, solution.p)
plant.animateTrajectory(timeStep, InitState, solution.tau, solution.x, solution.p)
solver.showResult(Info)
