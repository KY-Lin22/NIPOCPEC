%% Solving OCPEC with NIPOCPEC solver
clear all
clc
delete Biped.gif
delete Gen_InitialGuess.mat
[~, ~, ~] = rmdir('autoGen_CodeFiles', 's');

%% create an dynamics system 
timeStep = 0.01;
plant = Biped(timeStep);

% addpath(genpath('./autoGen_CodeFiles'))
plant.codeGen();

%% formulate init and reference state
nStages = 100;
% q = [x; z; torso; thigh_1; calf_1; thigh_2; calf_2] 
l_torso = plant.linkLength(1);
l_thigh = plant.linkLength(2);
l_calf = plant.linkLength(3);

% init state
x_Init = 0.2;
xita_torso_Init = -1/100*pi;
xita_thigh_1_Init = 1/25*pi;
xita_calf_1_Init = -1/75 * pi;
xita_thigh_2_Init = -1/40 * pi;

q_Init = plant.setInitConfiguration(x_Init, xita_torso_Init,...
    xita_thigh_1_Init, xita_calf_1_Init, xita_thigh_2_Init);
% q_Init(7) = -q_Init(7);
InitState = [q_Init;...
    0; 0; 0; 0; 0; 0; 0];

pf1_Init = q_Init(1) + l_thigh * sin(q_Init(4)) + l_calf * sin(q_Init(5));
pf2_Init = q_Init(1) + l_thigh * sin(q_Init(6)) + l_calf * sin(q_Init(7));
strd = 2 * abs(pf1_Init - pf2_Init);

% mid state
x_Mid = x_Init + 1/2 * strd;
xita_torso_Mid = 0;
xita_thigh_1_Mid = q_Init(6);
xita_calf_1_Mid = q_Init(7);
xita_thigh_2_Mid = q_Init(4);

q_Mid = plant.setInitConfiguration(x_Mid, xita_torso_Mid,...
    xita_thigh_1_Mid, xita_calf_1_Mid, xita_thigh_2_Mid);
% q_Mid(7) = -q_Mid(7);
MidState = [q_Mid;...
    0; 0; 0; 0; 0; 0; 0];

% final state
q_End = q_Init;
q_End(1) = q_End(1) + strd;

EndState = [q_End;...
    0; 0; 0; 0; 0; 0; 0];

pf1_End = q_End(1) + l_thigh * sin(q_End(4)) + l_calf * sin(q_End(5));
pf2_End = q_End(1) + l_thigh * sin(q_End(6)) + l_calf * sin(q_End(7));

% show initial and reference configuration of given plant
plant.plotConfiguration(InitState, [InitState, MidState, EndState])

%
tau_Max = [30; 30; 30; 30; 30;...
    1000; 1000];
tau_Min = [-30; -30; -30; -30; -30;...
    -1000; -1000];
x_Max = [q_End(1) + 0.2; 0.6; q_Init(3) + 1/50*pi; q_Init(4) + 1/2 * pi; 0; q_Init(6) + 1/2 * pi; 0;...
    20; 20; 20; 20; 20; 20; 20];
x_Min = [0; 0.3; q_Init(3) - 1/50*pi; q_Init(4) - 1/6 * pi; q_Init(5) - 1/6 * pi; q_Init(6) - 1/6 * pi; q_Init(7) - 1/6 * pi;...
    -20; -20; -20; -20; -20; -20; -20];
plant.setDynVarLimit(tau_Max, tau_Min, x_Max, x_Min) ;

%% formulate an OCPEC problem 
xRef_init_mid = TrajectoryInterpolation(InitState, MidState, 50);
xRef_mid_end = TrajectoryInterpolation(MidState, EndState, 50);

StageCost.xRef = [xRef_init_mid, xRef_mid_end];
StageCost.tauRef = repmat([0; 0; 0; 0; 0; 0; 0], 1, nStages);
StageCost.xWeight = [70; 30; 30; 100; 100; 100; 100;...
    0.1; 0.1; 1; 0.1; 0.1; 0.1; 0.1];
StageCost.tauWeight = [0.1; 0.1; 0.1; 0.1; 0.1; 0.001; 0.001]; 

TerminalCost.xRef = EndState;
TerminalCost.tauRef = [0; 0; 0; 0; 0; 0; 0];
TerminalCost.xWeight = [30; 100; 50; 100; 100; 100; 100;...
    0.1; 0.1; 1; 0.1; 0.1; 0.1; 0.1];
TerminalCost.tauWeight = [0.1; 0.1; 0.1; 0.1; 0.1; 0.001; 0.001];

OCPEC = OCPEC_Perturbed(plant, timeStep, nStages, InitState, StageCost, TerminalCost, 'Reg_Scholtes');% 'SmoothingEquation', 'Reg_NCPs', 'Reg_Scholtes' 

% set C
W_T = [1, 0, 0, plant.linkLength(2) * cos(OCPEC.x(4)), plant.linkLength(3) * cos(OCPEC.x(5)), 0, 0;...
       1, 0, 0, 0, 0, plant.linkLength(2) * cos(OCPEC.x(6)), plant.linkLength(3) * cos(OCPEC.x(7))];
vel_T = W_T * OCPEC.x(8:end);
C_auxiVar = [OCPEC.p(2) - OCPEC.p(3) - vel_T(1);...
             OCPEC.p(5) - OCPEC.p(6) - vel_T(2)];
    
C = C_auxiVar;
OCPEC.setEqualityConstraints(C);

% set G
G_nonSlip = [0.01 - OCPEC.p(1)*vel_T(1);...
     0.01 - OCPEC.p(4)*vel_T(2)]; % penalty slip motion
 
foot_1_X = OCPEC.x(1) + l_thigh * sin(OCPEC.x(4)) + l_calf * sin(OCPEC.x(5));
foot_1_Y = OCPEC.x(2) - l_thigh * cos(OCPEC.x(4)) - l_calf * cos(OCPEC.x(5));
foot_2_X = OCPEC.x(1) + l_thigh * sin(OCPEC.x(6)) + l_calf * sin(OCPEC.x(7));
foot_2_Y = OCPEC.x(2) - l_thigh * cos(OCPEC.x(6)) - l_calf * cos(OCPEC.x(7));

a = 1/2 * strd;
b_upper = 0.040;
b_lower = 0.025;
focus_foot_1_X = pf1_Init + a;
focus_foot_2_X = pf2_Init + a;

foot_1_EllipseTraj_upper = (foot_1_X - focus_foot_1_X)^2/(a^2) + (foot_1_Y)^2/b_upper^2 - 1;
foot_2_EllipseTraj_upper = (foot_2_X - focus_foot_2_X)^2/(a^2) + (foot_2_Y)^2/b_upper^2 - 1;

foot_1_EllipseTraj_lower = (foot_1_X - focus_foot_1_X)^2/(a^2) + (foot_1_Y)^2/b_lower^2 - 1;
foot_2_EllipseTraj_lower = (foot_2_X - focus_foot_2_X)^2/(a^2) + (foot_2_Y)^2/b_lower^2 - 1;

G_foot_EllipseTraj_upper = -[foot_1_EllipseTraj_upper;...
                             foot_2_EllipseTraj_upper];  
G_foot_EllipseTraj_lower = [foot_1_EllipseTraj_lower;...
                            foot_2_EllipseTraj_lower];                  
                  
G = [G_nonSlip;...
    G_foot_EllipseTraj_upper;...
    G_foot_EllipseTraj_lower]; 
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
solver.Option.Tolerance.KKT_Error_Feasibility = 1e-1;
solver.Option.Tolerance.KKT_Error_Stationarity = 1e-1;

solver.Option.RegularParam.nu_J = 1e-7;
solver.Option.RegularParam.nu_G = 1e-7;
solver.Option.RegularParam.nu_H = 0;

solver.Option.LineSearch.stepSize_Min = 0.01;
solver.Option.employFeasibilityRestorationPhase = true;

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
robustTest_Num = 10;
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
