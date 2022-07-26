%% Solving OCPEC with NIPOCPEC solver
clear all
clc
delete Quadruped.gif
delete Gen_InitialGuess.mat
[~, ~, ~] = rmdir('autoGen_CodeFiles', 's');

%% create an dynamics system  
timeStep = 0.01;
plant = Quadruped(timeStep);
tau_Max = [33.5; 33.5; 33.5; 33.5; 33.5; 33.5; 33.5; 33.5;...
           1000; 1000; 1000; 1000];
tau_Min = [-33.5; -33.5; -33.5; -33.5; -33.5; -33.5; -33.5; -33.5;...
          -1000; -1000; -1000; -1000];
x_Max = [1; 1.2; 4*pi; 4*pi; 4*pi; 4*pi; 4*pi; 4*pi; 4*pi; 4*pi; 4*pi;...
         100; 100; 100; 100; 100; 100; 100; 100; 100; 100; 100];
x_Min = [0; 0; -4*pi; -4*pi; -4*pi; -4*pi; -4*pi; -4*pi; -4*pi; -4*pi; -4*pi;...
        -100; -100; -100; -100; -100; -100; -100; -100; -100; -100; -100];
plant.setDynVarLimit(tau_Max, tau_Min, x_Max, x_Min) ;
% addpath(genpath('./autoGen_CodeFiles'))
plant.codeGen();

%% formulate init and reference state
nStages = 100;
% q = [x; z; torso; thigh_1; calf_1; thigh_2; calf_2; thigh_3; calf_3; thigh_4; calf_4] 
% init state
xita_Init = 2/5 * pi;
x_Init = 0.6;
q_Init = plant.setInitConfiguration(x_Init, xita_Init);
InitState = [q_Init;...
    0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];

xita_Mid = 2/5 * pi;
x_Mid = x_Init;
q_Mid = plant.setInitConfiguration(x_Mid, xita_Mid);
q_Mid(1) = q_Mid(1) - 0.5 * plant.linkLength(1);
q_Mid(2) = q_Mid(2) + 0.5;
q_Mid(3) = q_Mid(3) + pi;
q_Mid(4) = q_Mid(4) + pi;
q_Mid(5) = q_Mid(5) + pi;
q_Mid(6) = q_Mid(6) + pi;
q_Mid(7) = q_Mid(7) + pi;
q_Mid(8) = q_Mid(8) + pi;
q_Mid(9) = q_Mid(9) + pi;
q_Mid(10) = q_Mid(10) + pi;
q_Mid(11) = q_Mid(11) + pi;
MidState = [q_Mid;...
    0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];

xita_Mid2 = 2/7 * pi;
x_Mid2 = x_Init;
q_Mid2 = plant.setInitConfiguration(x_Mid2, xita_Mid2);
q_Mid2(1) = q_Mid2(1) - plant.linkLength(1);
q_Mid2(3) = q_Mid2(3) + 2 * pi;
q_Mid2(4) = q_Mid2(4) + 2 * pi;
q_Mid2(5) = q_Mid2(5) + 2 * pi;
q_Mid2(6) = q_Mid2(6) + 2 * pi;
q_Mid2(7) = q_Mid2(7) + 2 * pi;
q_Mid2(8) = q_Mid2(8) + 2 * pi;
q_Mid2(9) = q_Mid2(9) + 2 * pi;
q_Mid2(10) = q_Mid2(10) + 2 * pi;
q_Mid2(11) = q_Mid2(11) + 2 * pi;
Mid2State = [q_Mid2;...
    0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];

xita_End = 2/7 * pi;
x_End = x_Init;
q_End = plant.setInitConfiguration(x_End, xita_End);
q_End(1) = q_End(1) - plant.linkLength(1);
q_End(3) = q_End(3) + 2 * pi;
q_End(4) = q_End(4) + 2 * pi;
q_End(5) = q_End(5) + 2 * pi;
q_End(6) = q_End(6) + 2 * pi;
q_End(7) = q_End(7) + 2 * pi;
q_End(8) = q_End(8) + 2 * pi;
q_End(9) = q_End(9) + 2 * pi;
q_End(10) = q_End(10) + 2 * pi;
q_End(11) = q_End(11) + 2 * pi;
EndState = [q_End;...
    0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];

% show initial and reference configuration of given plant
plant.plotConfiguration(InitState, [InitState, MidState, Mid2State, EndState])

%%  formulate an OCPEC problem 
xRef_init_mid = TrajectoryInterpolation(InitState, MidState, 40);
xRef_mid_mid2 = TrajectoryInterpolation(MidState, Mid2State, 40);
xRef_mid2_end = TrajectoryInterpolation(Mid2State, EndState, 20);

StageCost.xRef = [xRef_init_mid, xRef_mid_mid2, xRef_mid2_end];
StageCost.tauRef = repmat(zeros(12, 1), 1, nStages);
StageCost.xWeight = [10; 10; 50; 50; 50; 50; 50; 50; 50; 50; 50;...
    0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1];
StageCost.tauWeight = [0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1;...
    0.001; 0.001; 0.001; 0.001]; 

TerminalCost.xRef = EndState;
TerminalCost.tauRef = zeros(12, 1);
TerminalCost.xWeight = [10; 10; 50; 50; 50; 50; 50; 50; 50; 50; 50;...
    0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1];
TerminalCost.tauWeight = [0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1;...
    0.001; 0.001; 0.001; 0.001];

OCPEC = OCPEC_Perturbed(plant, timeStep, nStages, InitState, StageCost, TerminalCost, 'Reg_Scholtes');% 'SmoothingEquation', 'Reg_NCPs', 'Reg_Scholtes'  

% set C
l_torso = plant.linkLength(1);
l_thigh = plant.linkLength(2);
l_calf = plant.linkLength(3);
W_T = [1, 0, 0, l_thigh * cos(OCPEC.x(4)), l_calf * cos(OCPEC.x(5)), 0, 0, 0, 0, 0, 0;...
       1, 0, 0, 0, 0, l_thigh * cos(OCPEC.x(6)), l_calf * cos(OCPEC.x(7)), 0, 0, 0, 0;...
       1, 0, l_torso * cos(OCPEC.x(3)), 0, 0, 0, 0, l_thigh * cos(OCPEC.x(8)), l_calf * cos(OCPEC.x(9)), 0, 0;...
       1, 0, l_torso * cos(OCPEC.x(3)), 0, 0, 0, 0, 0, 0, l_thigh * cos(OCPEC.x(10)), l_calf * cos(OCPEC.x(11))];
vel_T = W_T * OCPEC.x(12:end);
C = [OCPEC.p(2) - OCPEC.p(3) - vel_T(1);...
    OCPEC.p(5) - OCPEC.p(6) - vel_T(2);...
    OCPEC.p(8) - OCPEC.p(9) - vel_T(3);...
    OCPEC.p(11) - OCPEC.p(12) - vel_T(4)];

OCPEC.setEqualityConstraints(C);

% set G
G_nonSlip = [0.001 - OCPEC.p(1)*vel_T(1);...
             0.001 - OCPEC.p(4)*vel_T(2);...
             0.001 - OCPEC.p(7)*vel_T(3);...
             0.001 - OCPEC.p(10)*vel_T(4)];
G_kine_torso_thigh = [OCPEC.x(2) - l_torso * cos(OCPEC.x(3));...% torso 
                      OCPEC.x(2) - l_thigh * cos(OCPEC.x(4));...% thigh 1
                      OCPEC.x(2) - l_thigh * cos(OCPEC.x(6));...% thigh 2
                      OCPEC.x(2) - l_torso * cos(OCPEC.x(3)) - l_thigh * cos(OCPEC.x(8));...% thigh 3
                      OCPEC.x(2) - l_torso * cos(OCPEC.x(3)) - l_thigh * cos(OCPEC.x(10))]; % thigh 4
G = [G_nonSlip;...
    G_kine_torso_thigh]; 

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

solver.Option.LineSearch.stepSize_Min = 0.001;
solver.Option.employFeasibilityRestorationPhase = true;

solver.Option.zInit = 1e-1; 
solver.Option.zEnd  = 1e-1;
solver.Option.sInit = 1e-1;
solver.Option.sEnd  = 1e-1;

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
