%% Solving OCPEC with Ipopt solver via CasADi
clc
clear all
delete Hopper.gif
[~, ~, ~] = rmdir('autoGen_CodeFiles', 's');
addpath('E:\GitHub\CasADi\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

timeStep = 0.01;
nStages = 200; 
s = 1e-7; % slack 
z = 1e-4;
%% Dynamics
% dynamics variables
tau_Dim = 3;% control dim
x_Dim = 8; % state dim
p_Dim = 3; % equilibrium var dim
w_Dim = p_Dim; % auxiliary var dim
u_Dim = tau_Dim + p_Dim + w_Dim; % u = [tau; p; w]

x = SX.sym('x', x_Dim);
u = SX.sym('u', u_Dim);
tau = u(1 : tau_Dim);
p = u(tau_Dim + 1 : tau_Dim + p_Dim);
w = u(tau_Dim + p_Dim + 1 : end);

% state equations
mb = 1;
ml = 0.1;
Ib = 0.25;
Il = 0.025;
mu = 1;
g = 9.8;

M = diag([mb + ml, mb + ml, Ib + Il, ml]);
C = [0;...
    (mb + ml)*g;...
    0;...
    0];% here C := C*dq - G
B = [0, -sin(x(3));...
    0, cos(x(3));...
    1, 0;...
    0, 1];
W_N = [0, 1, x(4)*sin(x(3)), -cos(x(3))];
W_T = [1, 0, x(4)*cos(x(3)), sin(x(3))];
H = -C + B * [tau(1); tau(2)] + W_N' * p(1) + W_T' * tau(3);
f = [x(5:8);...
    inv(M)*H];
f_func = Function('f',{x,u}, {f}, {'x', 'u'}, {'f'});

% equilibrium dynamics
EqlbmDyn_l = [0; 0; 0];
EqlbmDyn_u = [Inf;Inf;Inf];
Gap = x(2) - x(4) * cos(x(3));
pT_lb = tau(3) - (-mu * p(1));
pT_ub = (mu * p(1)) - tau(3);
K = [Gap;...
     pT_lb;...
     pT_ub];
K_func = Function('K',{x,u}, {K}, {'x', 'u'}, {'K'}); 

% reformulate equilibrium dynamics as a set of inequality and equality constriants using Scholtes reformulation
BVI = SX.sym('BVI', 4 * p_Dim, 1);
lbg_BVI = zeros(4 * p_Dim, 1);
ubg_BVI = zeros(4 * p_Dim, 1);
for i = 1 : p_Dim
    BVI(4 * (i - 1) + 1 : 4 * i) = [p(i);...
                                    w(i);...
                                    w(i) - K(i);...
                                    s - p(i) * w(i)];
    lbg_BVI(4 * (i - 1) + 1 : 4 * i) = [0; 0; 0; 0];
    ubg_BVI(4 * (i - 1) + 1 : 4 * i) = [Inf; Inf; 0; Inf];
end
BVI_func = Function('BVI', {x,u}, {BVI}, {'x', 'u'}, {'BVI'});

% equality constraint and inequality constraint
vel_T = W_T * x(5:8);
EqCstr = p(2) - p(3) - vel_T; % using two auxilary variable for vel_T to reformulate friction 
lbg_EqCstr = 0;
ubg_EqCstr = 0;
EqCstr_func = Function('EqCstr', {x,u}, {EqCstr}, {'x', 'u'}, {'EqCstr'});

IneqCstr = 0.01 - p(1) * vel_T; % penalty slip motion
lbg_IneqCstr = 0;
ubg_IneqCstr = Inf;
IneqCstr_func = Function('IneqCstr', {x,u}, {IneqCstr}, {'x', 'u'}, {'IneqCstr'});

%% OCPEC
% specify initial and end state, cost ref and weight matrix
InitState = [0.1; 0.5; 0; 0.5; 0; 0; 0; 0];
MidState1 = [0.3; 0.65; 0; 0.2; 0; 0; 0; 0];
MidState2 = [0.4; 0.5; 0; 0.5; 0; 0; 0; 0];
MidState3 = [0.6; 0.65; 0; 0.2; 0; 0; 0; 0];
RefState  = [0.7; 0.5; 0; 0.5; 0; 0; 0; 0];


xRef_init_mid1 = TrajectoryInterpolation(InitState, MidState1, 40);
xRef_mid1_mid2 = TrajectoryInterpolation(MidState1, MidState2, 40);
xRef_mid2_mid3 = TrajectoryInterpolation(MidState2, MidState3, 40);
xRef_mid3_end = TrajectoryInterpolation(MidState3, RefState, 40);
xRef_end_end = TrajectoryInterpolation(RefState, RefState, 40);

StageCost.xRef = [xRef_init_mid1, xRef_mid1_mid2, xRef_mid2_mid3, xRef_mid3_end, xRef_end_end];
StageCost.tauRef = repmat([0; 0; 0], 1, nStages);
StageCost.xWeight = [100; 100; 100; 100; 0.1; 0.1; 0.1; 0.1];
StageCost.tauWeight = [0.1; 0.1; 0.001]; 


TerminalCost.xRef = RefState;
TerminalCost.tauRef = [0; 0; 0];
TerminalCost.xWeight = [100; 100; 100; 100; 0.1; 0.1; 0.1; 0.1];
TerminalCost.tauWeight = [0; 0; 0];

pWeight = 0.001 * eye(p_Dim);
wWeight = 0.001 * eye(w_Dim);

% optimal variables and their bounds
X = SX.sym('X', x_Dim, nStages);
U = SX.sym('U', u_Dim, nStages); 

tau_Max = [50; 50; 100];
tau_Min = [-50; -50; -100];
x_Max = [0.8; 0.7; pi; 0.5; 10; 10; 5; 5];
x_Min = [0; 0; -pi; 0.2; -10; -10; -5; -5];

lbx = -Inf * ones(x_Dim + u_Dim, nStages);
ubx = Inf * ones(x_Dim + u_Dim, nStages);
lbx(1 : x_Dim + tau_Dim, :) = repmat([x_Min; tau_Min], 1, nStages);
ubx(1 : x_Dim + tau_Dim, :) = repmat([x_Max; tau_Max], 1, nStages);

% cost function and constraints
L = 0; % initial cost function
g_Dim = size([f; BVI; EqCstr; IneqCstr], 1);
g = SX.sym('g', g_Dim, nStages); % constraint function
lbg = zeros(g_Dim, nStages);
ubg = zeros(g_Dim, nStages);

lbg(size(f, 1) + 1 : end, :) = repmat([lbg_BVI; lbg_EqCstr; lbg_IneqCstr], 1, nStages);
ubg(size(f, 1) + 1 : end, :) = repmat([ubg_BVI; ubg_EqCstr; ubg_IneqCstr], 1, nStages);

for n = 1 : nStages
    if n == 1
        x_nPrev = InitState;
    else
        x_nPrev = X(:, n - 1);
    end
    x_n = X(:, n);
    u_n = U(:, n);
    tau_n = u_n(1 : tau_Dim, 1);
    p_n = u_n(tau_Dim + 1 : tau_Dim + p_Dim, 1);
    w_n = u_n(tau_Dim + p_Dim + 1 : end, 1);
    
    % cost function
    L_n = 0.5 * (x_n - StageCost.xRef(:, n))' * diag(StageCost.xWeight) * (x_n - StageCost.xRef(:, n))...
        + 0.5 * (tau_n - StageCost.tauRef(:, n))' * diag(StageCost.tauWeight) * (tau_n - StageCost.tauRef(:, n))...
        + 0.5 * p_n' * pWeight * p_n...
        + 0.5 * w_n' * wWeight * w_n;
    L = L + L_n * timeStep;
    if n == nStages
        L_Terminal = 0.5 * (x_n - TerminalCost.xRef)' * diag(TerminalCost.xWeight) * (x_n - TerminalCost.xRef)...
            + 0.5 * (tau_n - TerminalCost.tauRef)' * diag(TerminalCost.tauWeight) * (tau_n - TerminalCost.tauRef);
        L = L + L_Terminal;
    end
    
    % discretize dynamics by implicit euler method
    F_n = x_nPrev - x_n + timeStep * f_func(x_n, u_n);
    % reformulated equilibrium constraint, equality constraint and inequality constraint
    BVI_n = BVI_func(x_n, u_n);
    EqCstr_n = EqCstr_func(x_n, u_n);
    IneqCstr_n = IneqCstr_func(x_n, u_n);
    % constraint function
    g(:, n) = [F_n;...
        BVI_n;...
        EqCstr_n;...
        IneqCstr_n]; 
end

%% Solver
% option
Option = struct;
Option.ipopt.max_iter = 500;
Option.ipopt.tol = 1e-2;
Option.ipopt.mu_target = 0.5 * (z)^2;

% reshape optimal variable and constraint
XU = reshape([X;U], (x_Dim + u_Dim) * nStages, 1);
g = reshape(g, g_Dim * nStages, 1);
lbx = reshape(lbx, (x_Dim + u_Dim) * nStages, 1);
ubx = reshape(ubx, (x_Dim + u_Dim) * nStages, 1);
lbg = reshape(lbg, g_Dim * nStages, 1);
ubg = reshape(ubg, g_Dim * nStages, 1);

%
robustTest_Num = 50;
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
% initial guess
x_Init = InitState;
x_End = TerminalCost.xRef;
x_0 = TrajectoryInterpolation(x_Init, x_End, nStages);% x
tau_Init = randn(tau_Dim, 1);
tau_End = TerminalCost.tauRef;
tau_0 = TrajectoryInterpolation(tau_Init, tau_End, nStages); % tau

XU_0 = ones(x_Dim + u_Dim, nStages);
XU_0(1 : x_Dim + tau_Dim, :) = [x_0; tau_0];
XU_0 = reshape(XU_0, (x_Dim + u_Dim) * nStages, 1);

RobustTestRecord.InitialGuess{i, 1} = XU_0;
%
OCPEC = struct('f', L, 'x', XU, 'g', g);
solver = nlpsol('solver', 'ipopt', OCPEC, Option);

totalTimeStart = tic;
solution = solver('x0', XU_0, 'lbx', lbx, 'ubx', ubx,...
             'lbg', lbg, 'ubg', ubg);
RobustTestRecord.solution{i, 1} = solution;         
totalTime = toc(totalTimeStart);

if strcmp(solver.stats.return_status, 'Solve_Succeeded')
    successCase = successCase + 1;
    RobustTestRecord.iterNum(successCase, 1) = solver.stats.iter_count;
    RobustTestRecord.totalTime(successCase, 1) = totalTime;
    RobustTestRecord.cost(successCase, 1) = full(solution.f);
    % compute constraint violation
    XU_Opt = reshape(full(solution.x), (x_Dim + u_Dim), nStages);
    x_Opt = XU_Opt(1 : x_Dim, :);
    u_Opt = XU_Opt(x_Dim + 1 : end, :);   
    tau_Opt = u_Opt(1 : tau_Dim, :);
    p_Opt = u_Opt(tau_Dim + 1 : tau_Dim + p_Dim, :);
    cstr = reshape(full(solution.g), g_Dim, nStages);
    
    Eq = cstr([11, 15, 19, 21], :); % w - K = 0, C
    DynF = cstr(1 : x_Dim, :); % 1 : 8       
    InEq = [[x_Opt; tau_Opt] - repmat([x_Min; tau_Min], 1, nStages);...
        repmat([x_Max; tau_Max], 1, nStages) - [x_Opt; tau_Opt];...
        cstr(22, :)];    
    K_Value = zeros(p_Dim, nStages);
    for n = 1 : nStages
        K_n = K_func(x_Opt(:, n), u_Opt(:, n));
        K_Value(:, n) = full(K_n);
    end
    lpu = cstr([9, 10, 13, 14, 17, 18], :);
    G_residual = min([zeros((2 * (x_Dim + tau_Dim) + 1) * nStages, 1), reshape(InEq, [], 1)], [], 2);
    pK_residual = min([zeros(2 * p_Dim * nStages, 1), reshape(lpu, [], 1)], [], 2); % 9, 10, 13, 14, 17, 18   
    
    Complementarity_pK = zeros(p_Dim, nStages);
    for n = 1 : nStages
        for j = 1 : p_Dim
            l_Vio = max(0, EqlbmDyn_l(j) - p_Opt(j, n));
            K_l_VioScale = min(1, max(0, p_Opt(j, n) - EqlbmDyn_l(j)));
            l_ComlVio = max(l_Vio, K_l_VioScale * max(K_Value(j, n), 0));
            u_Vio = max(0, p_Opt(j, n) - EqlbmDyn_u(j));
            K_u_VioScale = min(1, max(0, EqlbmDyn_u(j) - p_Opt(j, n)));
            u_ComlVio = max(u_Vio, K_u_VioScale * max(-K_Value(j, n), 0));
            Complementarity_pK(j, n) = max(l_ComlVio, u_ComlVio);
        end
    end
    
    RobustTestRecord.eqCstr(successCase, 1) = max([norm(reshape(Eq, [], 1), Inf), norm(reshape(DynF, [], 1), Inf)]);
    RobustTestRecord.ineqCstr(successCase, 1) = max(norm(G_residual, Inf), norm(reshape(pK_residual, [], 1), Inf));
    RobustTestRecord.compCstr(successCase, 1)  = norm(reshape(Complementarity_pK, [], 1), Inf);   
end

end
save('RobustTest_IPOPT_Data.mat', 'RobustTestRecord');
disp('robustTest')
disp(['success/total: ', num2str(successCase), '/', num2str(robustTest_Num)])
disp(['time per iter: ', num2str(1000 * sum(RobustTestRecord.totalTime) /sum(RobustTestRecord.iterNum), '%10.3f'), ' ms/Iter' ])
disp(['iterations: ', num2str(sum(RobustTestRecord.iterNum) / successCase, '%10.3f'), '(mean); ',...
    num2str(max(RobustTestRecord.iterNum(1 : successCase, 1))), '(max); ',...
    num2str(min(RobustTestRecord.iterNum(1 : successCase, 1))), '(min)'])
disp(['totalTime [s]: ', num2str(sum(RobustTestRecord.totalTime) / successCase, '%10.3f'), '(mean); ',...
    num2str(max(RobustTestRecord.totalTime(1 : successCase, 1)), '%10.3f'), '(max); ',...
    num2str(min(RobustTestRecord.totalTime(1 : successCase, 1)), '%10.3f'), '(min)'])
disp(['cost: ', num2str(sum(RobustTestRecord.cost) / successCase, '%10.3f'), '(mean); ',...
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

%% plot simulation result
% XU_Opt = reshape(full(solution.x), (x_Dim + u_Dim), nStages);
% x_Opt = XU_Opt(1 : x_Dim, :);
% u_Opt = XU_Opt(x_Dim + 1 : end, :);
% optSolution.tau = u_Opt(1 : tau_Dim, :);
% optSolution.x = x_Opt;
% optSolution.p = u_Opt(tau_Dim + 1 : tau_Dim + p_Dim, :);
% 
% plant = Hopper(timeStep, [], [], []);
% plant.setDynVarLimit(tau_Max, tau_Min, x_Max, x_Min) ;
% plant.codeGen();
% plant.plotSimuResult(timeStep, InitState, optSolution.tau, optSolution.x, optSolution.p)
% plant.animateTrajectory(timeStep, InitState, optSolution.tau, optSolution.x, optSolution.p)