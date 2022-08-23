%% Solving OCPEC with NIPOCPEC solver
clear all
clc
delete Gen_InitialGuess.mat
[~, ~, ~] = rmdir('autoGen_CodeFiles', 's');

%% create an dynamics system
timeStep = 0.02;
plant = FilippovDI(timeStep);
plant.codeGen();

%% formulate an OCPEC problem 
nStages = 100;
InitState = -1;

StageCost.xRef = zeros(1, nStages);
StageCost.tauRef = zeros(1, nStages);
StageCost.xWeight = 2;
StageCost.tauWeight = 0.00001; 
TerminalCost.xRef = 5/3;
TerminalCost.tauRef = 0;
TerminalCost.xWeight = 2;
TerminalCost.tauWeight = 0.00001;

OCPEC = OCPEC_Perturbed(plant, timeStep, nStages, InitState, StageCost, TerminalCost, 'Reg_Scholtes');% 'SmoothingEquation', 'Reg_NCPs', 'Reg_Scholtes' 

% set C
C = OCPEC.x - (OCPEC.p(2) - OCPEC.p(1));
OCPEC.setEqualityConstraints(C);

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

solver.Option.RegularParam.nu_J = 1e-7;
solver.Option.RegularParam.nu_G = 1e-7;
solver.Option.RegularParam.nu_H = 0;

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

%% solving OCPEC (sensitivity test)
x0 = -1.9 : 0.02 : -0.9;
sEnd = [1e-1, 1e-2, 1e-3, 1e-4];
num_x2 = size(x0, 2);
num_sEnd = size(sEnd, 2);
costRecord = zeros(num_sEnd, num_x2);

Test_Num = 5;
solver.Option.printLevel = 0;
for k = 1 : num_sEnd
    solver.Option.sEnd = sEnd(k);
    for j = 1 : num_x2
        InitState = x0(j);
        solver.OCPEC.InitState = InitState;
        
        TestRecord.InitialGuess = cell(Test_Num, 1);
        TestRecord.solution = cell(Test_Num, 1);
        successCase = 0;
        TestRecord.iterNum = zeros(Test_Num, 1);
        TestRecord.totalTime = zeros(Test_Num, 1);
        TestRecord.cost = zeros(Test_Num, 1);
        
        for i = 1 : Test_Num
            solver.generateInitialGuess();
            Gen_InitialGuess = load('Gen_InitialGuess.mat');
            [solution, Info] = solver.solveOCPEC(Gen_InitialGuess.Iterate);
            TestRecord.InitialGuess{i, 1} = Gen_InitialGuess.Iterate;
            TestRecord.solution{i, 1} = solution;          
            if Info.iterProcess.terminalStatus == 1
                successCase = successCase + 1;
                TestRecord.iterNum(successCase, 1) = Info.iterProcess.iterNum;
                TestRecord.totalTime(successCase, 1) = Info.iterProcess.Time.total;
                TestRecord.cost(successCase, 1) = Info.solutionMsg.totalCost;
            end
            disp(['success / Test No.: ', num2str(successCase), ' / ', num2str(i)])
        end
        costRecord(k, j) = min(TestRecord.cost(1 : successCase, 1));
    end
        
end

figure(1)
plot(x0, costRecord(1, :), 'k',...
    x0, costRecord(2, :), 'b',...
    x0, costRecord(3, :), 'r',...
    x0, costRecord(4, :), 'g', 'LineWidth', 1.2)
legend('1e-1', '1e-2', '1e-3', '1e-4')
xlabel('x0')
ylabel('cost')
