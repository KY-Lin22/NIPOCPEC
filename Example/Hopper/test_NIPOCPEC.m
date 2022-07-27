%% Solving OCPEC with NIPOCPEC solver
clear all
clc
delete Hopper.gif
delete Gen_InitialGuess.mat
[~, ~, ~] = rmdir('autoGen_CodeFiles', 's');

%% create an dynamics system 
timeStep = 0.01;
mass = [1; 0.1];
inertia = [0.25; 0.025];
mu = 0.45;
plant = Hopper(timeStep, mass, inertia, mu);
tau_Max = [50; 50; 100];
tau_Min = [-50; -50; -100];
x_Max = [0.8; 1.5; pi; 0.5; 10; 10; 5; 5];
x_Min = [0; 0; -pi; 0.1; -10; -10; -5; -5];
plant.setDynVarLimit(tau_Max, tau_Min, x_Max, x_Min) ;
plant.codeGen();

%% formulate an OCPEC problem 
nStages = 100;
InitState = [0.1; 0.5; 0; 0.5; 0; 0; 0; 0];
MidState = [0.4; 0.8; 0; 0.1; 0; 0; 0; 0];
RefState = [0.7; 0.5; 0; 0.5; 0; 0; 0; 0];

xRef_init_mid = TrajectoryInterpolation(InitState, MidState, 50);
xRef_mid_end = TrajectoryInterpolation(MidState, RefState, 50);
StageCost.xRef = [xRef_init_mid, xRef_mid_end];
StageCost.tauRef = repmat([0; 0; 0], 1, nStages);
StageCost.xWeight = [50; 50; 20; 50; 0.1; 0.1; 0.1; 0.1];
StageCost.tauWeight = [0.1; 0.1; 0.001]; 

TerminalCost.xRef = RefState;
TerminalCost.tauRef = [0; 0; 0];
TerminalCost.xWeight = [50; 50; 50; 50; 0.1; 0.1; 0.1; 0.1];
TerminalCost.tauWeight = [0.1; 0.1; 0.001];

% show initial and reference configuration of given plant
plotConfiguration(plant, InitState, [InitState, MidState, RefState])

OCPEC = OCPEC_Perturbed(plant, timeStep, nStages, InitState, StageCost, TerminalCost, 'Reg_Scholtes');% 'SmoothingEquation', 'Reg_NCPs', 'Reg_Scholtes'  

W_T = [1, 0, OCPEC.x(4)*cos(OCPEC.x(3)), sin(OCPEC.x(3))];
vel_T = W_T * OCPEC.x(5:8);
C = OCPEC.p(2) - OCPEC.p(3) - vel_T; % using two auxilary variable for vel_T to reformulate friction 

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

%% solving OCPEC (robust test)
robustTest_Num = 100;
solver.Option.printLevel = 0;
RobustTestRecord.InitialGuess = cell(robustTest_Num, 1);
RobustTestRecord.solution = cell(robustTest_Num, 1);
successCase = 0;
RobustTestRecord.iterNum = zeros(robustTest_Num, 1);
RobustTestRecord.totalTime = zeros(robustTest_Num, 1);
RobustTestRecord.cost = zeros(robustTest_Num, 1);
RobustTestRecord.eqCstr = zeros(robustTest_Num, 1);
RobustTestRecord.ineqCstr = zeros(robustTest_Num, 1);
RobustTestRecord.compCstr = zeros(robustTest_Num, 1);

for i = 1 : robustTest_Num    
    
    solver.generateInitialGuess();
    Gen_InitialGuess = load('Gen_InitialGuess.mat');    
    [solution, Info] = solver.solveOCPEC(Gen_InitialGuess.Iterate);
    RobustTestRecord.InitialGuess{i, 1} = Gen_InitialGuess.Iterate;
    RobustTestRecord.solution{i, 1} = solution;

    if Info.iterProcess.terminalStatus == 1   
        successCase = successCase + 1;
        RobustTestRecord.iterNum(successCase, 1) = Info.iterProcess.iterNum;
        RobustTestRecord.totalTime(successCase, 1) = Info.iterProcess.Time.total;
        
        RobustTestRecord.cost(successCase, 1) = Info.solutionMsg.totalCost;    
        
        RobustTestRecord.eqCstr(successCase, 1) = max([Info.solutionMsg.r_eq_C, Info.solutionMsg.r_eq_F]);
        RobustTestRecord.ineqCstr(successCase, 1) = max([Info.solutionMsg.r_ineq_G, Info.solutionMsg.r_eqlb_ineq]);
        RobustTestRecord.compCstr(successCase, 1)  = Info.solutionMsg.r_eqlb_comp;
    end
    disp(['success / Test No.: ', num2str(successCase), ' / ', num2str(i)])
end  
% show result 
disp('robustTest')
disp(['success/total: ', num2str(successCase), '/', num2str(robustTest_Num)])
disp(['time per iter: ', num2str(1000 * sum(RobustTestRecord.totalTime) /sum(RobustTestRecord.iterNum), '%10.3f'), ' ms/Iter' ])
disp(['iterations: ', num2str(sum(RobustTestRecord.iterNum) / successCase, '%10.3f'), '(mean); ',...
    num2str(max(RobustTestRecord.iterNum(1 : successCase, 1))), '(max); ',...
    num2str(min(RobustTestRecord.iterNum(1 : successCase, 1))), '(min)'])
disp(['totalTime [s]: ', num2str(sum(RobustTestRecord.totalTime) / successCase, '%10.3f'), '(mean); ',...
    num2str(max(RobustTestRecord.totalTime(1 : successCase, 1)), '%10.3f'), '(max); ',...
    num2str(min(RobustTestRecord.totalTime(1 : successCase, 1)), '%10.3f'), '(min)'])
disp(['cost: ', num2str(sum(RobustTestRecord.cost) / successCase , '%10.3f'), '(mean); ',...
    num2str(max(RobustTestRecord.cost(1 : successCase, 1)), '%10.3f'), '(max); ',...
    num2str(min(RobustTestRecord.cost(1 : successCase, 1)), '%10.3f'),'(min)'])
disp(['eqCstr: ', num2str(sum(RobustTestRecord.eqCstr) / successCase,'%10.3e'), '(mean); ',...
    num2str(max(RobustTestRecord.eqCstr(1 : successCase, 1)), '%10.3e'), '(max); ',...
    num2str(min(RobustTestRecord.eqCstr(1 : successCase, 1)), '%10.3e'),'(min)'])
disp(['ineqCstr: ', num2str(sum(RobustTestRecord.ineqCstr) / successCase, '%10.3e'), '(mean); ',...
    num2str(max(RobustTestRecord.ineqCstr(1 : successCase, 1)), '%10.3e'), '(max); ',...
    num2str(min(RobustTestRecord.ineqCstr(1 : successCase, 1)), '%10.3e'),'(min)'])
disp(['compCstr: ', num2str(sum(RobustTestRecord.compCstr) / successCase, '%10.3e'), '(mean); ',...
    num2str(max(RobustTestRecord.compCstr(1 : successCase, 1)), '%10.3e'), '(max); ',...
    num2str(min(RobustTestRecord.compCstr(1 : successCase, 1)), '%10.3e'),'(min)'])
