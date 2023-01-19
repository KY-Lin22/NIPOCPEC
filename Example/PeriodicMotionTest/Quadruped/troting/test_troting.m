%% Solving OCPEC with NIPOCPEC solver
clear all
clc
delete Quadruped.gif
delete Gen_InitialGuess.mat
[~, ~, ~] = rmdir('autoGen_CodeFiles', 's');

%% create an dynamics system  
PeriodicDVI_InitGuess = load('PeriodicDVI_InitGuess.mat');
timeStep = 0.01;
plant = Quadruped_PeriodicGait(timeStep);
plant.periodLimitCycle = PeriodicDVI_InitGuess.OscillatorParam.periodLimitCycle;
plant.dutyFactor = PeriodicDVI_InitGuess.OscillatorParam.dutyFactor;
plant.InitPhase = PeriodicDVI_InitGuess.OscillatorParam.InitPhase;
plant.k_N = PeriodicDVI_InitGuess.FootVelFuncParam.k_N;
plant.k_T = PeriodicDVI_InitGuess.FootVelFuncParam.k_T;
%%
plant.codeGen();
%% formulate init and reference state
l_torso = plant.linkLength(1);
l_thigh = plant.linkLength(2);
l_calf = plant.linkLength(3);

nStages = 25;
% q = [x; z; torso; thigh_1; calf_1; thigh_2; calf_2; thigh_3; calf_3; thigh_4; calf_4] 
% init state
xita_Init = 1/5 * pi;
x_Init = PeriodicDVI_InitGuess.x_Init;
q_Init = plant.setInitConfiguration(x_Init, xita_Init);
x_PrdDVI_Init = PeriodicDVI_InitGuess.InitState;
InitState = [q_Init; zeros(plant.qDim, 1); x_PrdDVI_Init];

% middle
desiredAng = pi/10;
q_Mid = [PeriodicDVI_InitGuess.solution.x(12, 25) - l_calf*sin(desiredAng) + l_thigh*sin(acos((q_Init(2) - l_calf * cos(desiredAng))/l_thigh));...% x: leg1_T - l3sinx5 + l2sinx4
    q_Init(2);...% y = q_Init's y
    pi/2;... % torso
    -acos((q_Init(2) - l_calf * cos(desiredAng))/l_thigh);...
    desiredAng;...
    -acos((q_Init(2) - l_calf * cos(xita_Init))/l_thigh);...
    xita_Init;...
    -acos((q_Init(2) - l_calf * cos(xita_Init))/l_thigh);...
    xita_Init;...
    -acos((q_Init(2) - l_calf * cos(desiredAng))/l_thigh);...
    desiredAng]; 
x_PrdDVI_MId = PeriodicDVI_InitGuess.solution.x(:, 25);
MidState = [q_Mid; zeros(plant.qDim, 1); x_PrdDVI_MId];

% end
xita_End = xita_Init;
x_End = PeriodicDVI_InitGuess.solution.x(12, nStages); 
q_End = plant.setInitConfiguration(x_End, xita_End);
x_PrdDVI_End = PeriodicDVI_InitGuess.solution.x(:, nStages);
EndState = [q_Mid; zeros(plant.qDim, 1); x_PrdDVI_End];

% show initial and reference configuration of given plant
plant.plotConfiguration(InitState, MidState)

%%  formulate an OCPEC problem 
% cost function
xRef_init_mid = TrajectoryInterpolation(InitState, MidState, 25);
% xRef_mid_end = TrajectoryInterpolation(MidState, EndState, 25);
xRef_init_end = xRef_init_mid;

xRef_init_end(23 : 38, :) = PeriodicDVI_InitGuess.solution.x(:, 1 : nStages);
%%
StageCost.xRef = xRef_init_end;
StageCost.tauRef = repmat(zeros(12, 1), 1, nStages);
StageCost.xWeight = [1; 1; 1; 0.5 * ones(8, 1);...% rigid body position
    0.1 * ones(11, 1);...% rigid body velocity
    100 * ones(16, 1)]; % periodic DVI
StageCost.tauWeight = [0.1 * ones(8, 1);...% control
    0.001 * ones(4, 1)]; % friction

TerminalCost.xRef = EndState;
TerminalCost.tauRef = zeros(12, 1);
TerminalCost.xWeight = [1; 1; 0.5; 0.5 * ones(8, 1);...% rigid body position
    0.1 * ones(11, 1);...% rigid body velocity
    100 * ones(16, 1)]; % periodic DVI
TerminalCost.tauWeight = [0.1 * ones(8, 1);...% control
    0.001 * ones(4, 1)]; % friction
OCPEC = OCPEC_Perturbed(plant, timeStep, nStages, InitState, StageCost, TerminalCost, 'Reg_Scholtes');% 'SmoothingEquation', 'Reg_NCPs', 'Reg_Scholtes'  

% set C
W_T = [1, 0, 0, l_thigh * cos(OCPEC.x(4)), l_calf * cos(OCPEC.x(5)), 0, 0, 0, 0, 0, 0;...
       1, 0, 0, 0, 0, l_thigh * cos(OCPEC.x(6)), l_calf * cos(OCPEC.x(7)), 0, 0, 0, 0;...
       1, 0, l_torso * cos(OCPEC.x(3)), 0, 0, 0, 0, l_thigh * cos(OCPEC.x(8)), l_calf * cos(OCPEC.x(9)), 0, 0;...
       1, 0, l_torso * cos(OCPEC.x(3)), 0, 0, 0, 0, 0, 0, l_thigh * cos(OCPEC.x(10)), l_calf * cos(OCPEC.x(11))];
vel_T = W_T * OCPEC.x(12:22);
C_auxiVar_vel_T = [OCPEC.p(2) - OCPEC.p(3) - vel_T(1);...
    OCPEC.p(5) - OCPEC.p(6) - vel_T(2);...
    OCPEC.p(8) - OCPEC.p(9) - vel_T(3);...
    OCPEC.p(11) - OCPEC.p(12) - vel_T(4)]; % auxiliary variable for tangential velocity

footPosiN_1 = OCPEC.x(2) - l_thigh * cos(OCPEC.x(4)) - l_calf * cos(OCPEC.x(5));
footPosiT_1 = OCPEC.x(1) + l_thigh * sin(OCPEC.x(4)) + l_calf * sin(OCPEC.x(5));
footPosiN_2 = OCPEC.x(2) - l_thigh * cos(OCPEC.x(6)) - l_calf * cos(OCPEC.x(7));
footPosiT_2 = OCPEC.x(1) + l_thigh * sin(OCPEC.x(6)) + l_calf * sin(OCPEC.x(7));
footPosiN_3 = OCPEC.x(2) - l_torso * cos(OCPEC.x(3)) - l_thigh * cos(OCPEC.x(8)) - l_calf * cos(OCPEC.x(9));
footPosiT_3 = OCPEC.x(1) + l_torso * sin(OCPEC.x(3)) + l_thigh * sin(OCPEC.x(8)) + l_calf * sin(OCPEC.x(9));
footPosiN_4 = OCPEC.x(2) - l_torso * cos(OCPEC.x(3)) - l_thigh * cos(OCPEC.x(10)) - l_calf * cos(OCPEC.x(11));
footPosiT_4 = OCPEC.x(1) + l_torso * sin(OCPEC.x(3)) + l_thigh * sin(OCPEC.x(10)) + l_calf * sin(OCPEC.x(11));

C_footPos = [OCPEC.x(25) - footPosiN_3; OCPEC.x(26) - footPosiT_3;...
    OCPEC.x(29) - footPosiN_4; OCPEC.x(30) - footPosiT_4;...
    OCPEC.x(33) - footPosiN_1; OCPEC.x(34) - footPosiT_1;...
    OCPEC.x(37) - footPosiN_2; OCPEC.x(38) - footPosiT_2];

C = [C_auxiVar_vel_T; C_footPos];
OCPEC.setEqualityConstraints(C);

% set G
safeHeight = 0.01;
G_kine_torso_thigh = [OCPEC.x(2) - l_torso * cos(OCPEC.x(3)) - safeHeight;...% torso 
                      OCPEC.x(2) - l_thigh * cos(OCPEC.x(4)) - safeHeight;...% thigh 1
                      OCPEC.x(2) - l_thigh * cos(OCPEC.x(6)) - safeHeight;...% thigh 2
                      OCPEC.x(2) - l_torso * cos(OCPEC.x(3)) - l_thigh * cos(OCPEC.x(8)) - safeHeight;...% thigh 3
                      OCPEC.x(2) - l_torso * cos(OCPEC.x(3)) - l_thigh * cos(OCPEC.x(10)) - safeHeight]; % thigh 4

G = G_kine_torso_thigh;           
OCPEC.setInequalityConstraints(G);

% generate necessary code files
OCPEC.codeGen()

OCPEC.showInfo()

%% create a NIPOCPEC_Solver object 
solver = NIPOCPEC_Solver(OCPEC);
% generate necessary code files
solver.codeGen();

%% set option and generate initial guess
solver.Option.maxIterNum = 800;
solver.Option.Tolerance.KKT_Error_Total = 1e-1;
solver.Option.Tolerance.KKT_Error_Feasibility = 1e-1;
solver.Option.Tolerance.KKT_Error_Stationarity = 1e-1;

solver.Option.RegularParam.nu_J = 1e-7;
solver.Option.RegularParam.nu_G = 1e-6;
solver.Option.RegularParam.nu_H = 0;

solver.Option.LineSearch.stepSize_Min = 0.001;
solver.Option.employFeasibilityRestorationPhase = true;
solver.Option.FRP.maxIterNum = 50;
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
IterateInit.x = StageCost.xRef;
[solution, Info] = solver.solveOCPEC(IterateInit);

%%
plant.animateTrajectory(timeStep, InitState(1: 22), solution.tau, solution.x(1: 22, :), solution.p(1 : 12, :))
