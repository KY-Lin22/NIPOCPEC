%% Solving OCPEC with Ipopt solver via CasADi
clc
clear all
delete AffineDVI.gif

addpath('E:\GitHub\CasADi\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

timeStep = 0.01; % Time horizon
nStages = 100; % number of control intervals

%% Dynamics
% dynamics variables
tau_Dim = 1;% control dim
x_Dim = 2; % state dim
p_Dim = 1; % equilibrium var dim
w_Dim = p_Dim; % auxiliary var dim
u_Dim = tau_Dim + p_Dim + w_Dim; % u = [tau; p; w]

x = SX.sym('x', x_Dim);
u = SX.sym('u', u_Dim);
tau = u(1 : tau_Dim);
p = u(tau_Dim + 1 : tau_Dim + p_Dim);
w = u(tau_Dim + p_Dim + 1 : end);

% system dynamics equations
A = [1, -3; ...
    -8, 10];
B = [-3;...
    -1];
F = [4;...
    8];
f = A * x + B * p + F * tau; % xDot = f(tau, x, p)
f_func = Function('f',{x,u}, {f}, {'x', 'u'}, {'f'});

% equilibrium dynamics equation
EqlbmDyn_l = -1;
EqlbmDyn_u = 1;
C = [1, -3];
D = 5;
E = 3;
K = C * x + D * p + E * tau;
K_func = Function('K',{x,u}, {K}, {'x', 'u'}, {'K'});

% reformulate equilibrium dynamics as a set of inequality and equality constriants using Scholtes reformulation
s = 1e-3; % slack 
BVI = [p - EqlbmDyn_l;...
    EqlbmDyn_u - p;...% G_BVI >= 0
    w - K;...% C_BVI = 0
    s - (p - EqlbmDyn_l) * w;...
    s + (EqlbmDyn_u - p) * w];% PHI >= 0
BVI_func = Function('BVI', {x,u}, {BVI}, {'x', 'u'}, {'BVI'});
lbg_BVI = [0; 0; 0; 0; 0];
ubg_BVI = [Inf; Inf; 0; Inf; Inf];

%% OCPEC
% specify initial and end state, cost ref and weight matrix
InitState = [-1/2; -1];
EndState = [0; 0];

StageCost.xRef = repmat(EndState, 1, nStages);
StageCost.tauRef = zeros(1, nStages);
StageCost.xWeight = [20; 20];
StageCost.tauWeight = 1;
TerminalCost.xRef = EndState;
TerminalCost.tauRef = 0;
TerminalCost.xWeight = [20; 20];
TerminalCost.tauWeight = 1;

pWeight = 0.001 * eye(p_Dim);
wWeight = 0.001 * eye(w_Dim);

% optimal variables and their bounds
X = SX.sym('X', x_Dim, nStages);
U = SX.sym('U', u_Dim, nStages); 

tau_Max = 2;
tau_Min = -2;
x_Max = [2; 2];
x_Min = [-2; -2];

lbx = -Inf * ones(x_Dim + u_Dim, nStages);
ubx = Inf * ones(x_Dim + u_Dim, nStages);
lbx(1 : x_Dim + tau_Dim, :) = repmat([x_Min; tau_Min], 1, nStages);
ubx(1 : x_Dim + tau_Dim, :) = repmat([x_Max; tau_Max], 1, nStages);

% cost function and constraints
L = 0; % initial cost function
g_Dim = size([f; BVI], 1);
g = SX.sym('g', g_Dim, nStages); % constraint function
lbg = zeros(g_Dim, nStages);
ubg = zeros(g_Dim, nStages);
lbg(size(f, 1) + 1 : end, :) = repmat(lbg_BVI, 1, nStages);
ubg(size(f, 1) + 1 : end, :) = repmat(ubg_BVI, 1, nStages);

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
    % reformulated equilibrium constraint
    BVI_n = BVI_func(x_n, u_n);
    
    % constraint function
    g(:, n) = [F_n;...
        BVI_n]; 
end

%% Solver
% option
Option = struct;
Option.ipopt.max_iter = 500;
Option.ipopt.tol = 1e-2;
z = 1e-3;
Option.ipopt.mu_target = 0.5 * (z)^2;

% reshape optimal variable and constraint
XU = reshape([X;U], (x_Dim + u_Dim) * nStages, 1);
g = reshape(g, g_Dim * nStages, 1);
lbx = reshape(lbx, (x_Dim + u_Dim) * nStages, 1);
ubx = reshape(ubx, (x_Dim + u_Dim) * nStages, 1);
lbg = reshape(lbg, g_Dim * nStages, 1);
ubg = reshape(ubg, g_Dim * nStages, 1);

%
robustTest_Num = 1;
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

XU_0 = zeros(x_Dim + u_Dim, nStages);
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
    RobustTestRecord.cost(i, 1) = full(solution.f);
    % compute constraint violation
    XU_Opt = reshape(full(solution.x), (x_Dim + u_Dim), nStages);
    x_Opt = XU_Opt(1 : x_Dim, :);
    u_Opt = XU_Opt(x_Dim + 1 : end, :);   
    tau_Opt = u_Opt(1 : tau_Dim, :);
    p_Opt = u_Opt(tau_Dim + 1 : tau_Dim + p_Dim, :);
    cstr = reshape(full(solution.g), g_Dim, nStages);
    EQ = cstr(5, :); % w - K = 0
    DynF = cstr(1 : x_Dim, :); % 1 : 2       
    G = [[x_Opt; tau_Opt] - repmat([x_Min; tau_Min], 1, nStages);...
        repmat([x_Max; tau_Max], 1, nStages) - [x_Opt; tau_Opt]];    
    K_Value = zeros(p_Dim, nStages);
    for n = 1 : nStages
        K_n = K_func(x_Opt(:, n), u_Opt(:, n));
        K_Value(:, n) = full(K_n);
    end
    lpu = cstr(x_Dim + 1: x_Dim + 2 * p_Dim, :);
    G_residual = min([zeros(2 * (x_Dim + tau_Dim) * nStages, 1), reshape(G, [], 1)], [], 2);
    pK_residual = min([zeros(2 * p_Dim * nStages, 1), reshape(lpu, [], 1)], [], 2); % 3 : 4   
    
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
    
    RobustTestRecord.eqCstr(successCase, 1) = max([norm(reshape(EQ, [], 1), Inf), norm(reshape(DynF, [], 1), Inf)]);
    RobustTestRecord.ineqCstr(successCase, 1) = max(norm(G_residual, Inf), norm(reshape(pK_residual, [], 1), Inf));
    RobustTestRecord.compCstr(successCase, 1)  = norm(reshape(Complementarity_pK, [], 1), Inf);   
end

end

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

% %% plot simulation result
% % extract solution
% XU_Opt = reshape(full(solution.x), (x_Dim + u_Dim), nStages);
% x_Opt = XU_Opt(1 : x_Dim, :);
% u_Opt = XU_Opt(x_Dim + 1 : end, :);
% K_Opt = zeros(p_Dim, nStages);
% for n = 1 : nStages
%     K_Opt_n = K_func(x_Opt(:, n), u_Opt(:, n));
%     K_Opt(:, n) = full(K_Opt_n);
% end
% 
% tau_Opt = u_Opt(1 : tau_Dim, :);
% p_Opt = u_Opt(tau_Dim + 1 : tau_Dim + p_Dim, :);
% 
% timeAxis = 0 : timeStep : nStages * timeStep;
% 
% figure(111)
% subplot(3,1,1)
% plot(timeAxis, [InitState(1), x_Opt(1, :)], 'r',...
%      timeAxis, [InitState(2), x_Opt(2, :)], 'g', 'LineWidth',1.2)
% legend('x1', 'x2')
% xlabel('time(s)')
% title('system state')
% 
% subplot(3,1,2)
% plot(timeAxis(2:end), tau_Opt(1,:), 'LineWidth', 1.2)
% xlabel('time(s)')
% title('control')
% 
% subplot(3,1,3)
% plot(timeAxis(2:end), p_Opt(1, :), 'k',...
%      timeAxis(2:end), K_Opt(1, :), 'b', 'LineWidth', 1.2)
% legend('p', 'K') 
% xlabel('time(s)')
% title('equilibrium dynamics')
% 
% %% animation trajectory
% % get baseline coordinate
% baseLine_horizontal_X = [x_Min(1); x_Max(1)];
% baseLine_horizontal_Y = [0; 0];  
% baseLine_vertical_X = [0; 0];
% baseLine_vertical_Y = [x_Min(2); x_Max(2)]; 
% 
% % get state sequence based on given x sequence
% state_X = [InitState(1), x_Opt(1, :)];
% state_Y = [InitState(2), x_Opt(2, :)];
% 
% % generate a trajectory-time sequence for animation
% traj_X = cell(1, nStages + 1);
% traj_Y = cell(1, nStages + 1);
% for n = 1 : nStages + 1
%     traj_X{1, n} = [state_X(1, 1 : n), repmat(state_X(1, n), 1, nStages + 1 - n)];
%     traj_Y{1, n} = [state_Y(1, 1 : n), repmat(state_Y(1, n), 1, nStages + 1 - n)];
% end
% 
% % compute value of p, K and generate a trajectory-time sequence for animation
% p_Opt = [p_Opt(:, 1), p_Opt];
% K_Opt = [K_Opt(:, 1), K_Opt];% extend
% trajp_X = cell(1, nStages + 1);
% trajp_Y = cell(1, nStages + 1);
% trajK_X = cell(1, nStages + 1);
% trajK_Y = cell(1, nStages + 1);
% for n = 1 : nStages + 1
%     trajp_X{1, n} = [timeAxis(1:n), repmat(timeAxis(n), 1, nStages + 1 - n)];
%     trajp_Y{1, n} = [p_Opt(1, 1:n), repmat(p_Opt(1, n), 1, nStages + 1 - n)];
%     trajK_X{1, n} = [timeAxis(1:n), repmat(timeAxis(n), 1, nStages + 1 - n)];
%     trajK_Y{1, n} = [K_Opt(1, 1:n), repmat(K_Opt(1, n), 1, nStages + 1 - n)];
% end
% 
% % get figure size
% figure(100)
% pos = get(gcf, 'Position');
% width = pos(3);
% height = pos(4);
% % define movie record
% mov = zeros(height*3/2, width*3/2, 1, nStages + 1, 'uint8');
% % pre allocate
% subplot(1,1,1)
% 
% plot(baseLine_horizontal_X, baseLine_horizontal_Y, '.-k', 'MarkerSize', 0.5, 'LineWidth', 0.5);
% hold on
% plot(baseLine_vertical_X, baseLine_vertical_Y, '.-k', 'MarkerSize', 0.5, 'LineWidth', 0.5);
% hold on
% state = plot(state_X(1, 1), state_Y(1, 1), '*', 'Color', [1 0 0], 'MarkerSize', 6, 'LineWidth', 2);
% hold on
% traj = plot(traj_X{1, 1}, traj_Y{1, 1}, '-', 'Color', [1 0 0], 'MarkerSize', 1, 'LineWidth', 1);
% 
% axisLimit_X = baseLine_horizontal_X;
% axisLimit_Y = baseLine_vertical_Y; 
% axis([axisLimit_X; axisLimit_Y]);
% % timePrint = title(sprintf('Time: %0.2f sec', timeAxis(1)));
% xlabel('x1')
% ylabel('x2')
% % subplot(2,1,2)
% % trajp = plot(trajp_X{1, 1}, trajp_Y{1, 1}, 'k', 'LineWidth', 1.2);
% % hold on
% % trajK = plot(trajK_X{1, 1}, trajK_Y{1, 1}, 'b', 'LineWidth', 1.2);
% % axis([0; timeAxis(end); min([p_Opt,K_Opt]); max([p_Opt,K_Opt])]);
% % legend('p', 'K')
% % xlabel('time(s)')
% % title('equilibrium dynamics');
% % animation
% for n = 1 : nStages + 1
%     % update XData and YData   
%     set(state, 'XData', state_X(1, n), 'YData', state_Y(1, n));
%     set(traj, 'XData', traj_X{1, n}, 'YData', traj_Y{1, n});
% %     set(timePrint, 'String', sprintf('Time: %0.2f sec', timeAxis(n)));
% %     set(trajp, 'XData', trajp_X{1, n}, 'YData', trajp_Y{1, n});
% %     set(trajK, 'XData', trajK_X{1, n}, 'YData', trajK_Y{1, n});  
%     % get frame as an image
%     f = getframe(gcf);
%     % Create a colormap for the first frame. for the rest of the frames, use the same colormap
%     if n == 1
%         [mov(:,:,1,n), map] = rgb2ind(f.cdata, 256, 'nodither');
%     else
%         mov(:,:,1,n) = rgb2ind(f.cdata, map, 'nodither');
%     end
% end
% % create an animated GIF
% imwrite(mov, map, 'AffineDVI.gif', 'DelayTime', 0, 'LoopCount', inf)
