function setDynamics(plant)
%setDynamics
%   Detailed explanation goes here

% rigid body dynamics
q = plant.x(1 : 11);
dq = plant.x(12 : 22);

tau_control = plant.tau(1 : 8);% control

pN_1 = plant.p(1);% normal contact force in leg 1
pT_1 = plant.tau(9); % friction force in leg 1

pN_2 = plant.p(4); % normal contact force in leg 2
pT_2 = plant.tau(10); % friction force in leg 2

pN_3 = plant.p(7);% normal contact force in leg 3
pT_3 = plant.tau(11); % friction force in leg 3

pN_4 = plant.p(10); % normal contact force in leg 4
pT_4 = plant.tau(12); % friction force in leg 4

mu = plant.mu;

% periodic DVI
x_PrdDVI = plant.x(23 : end);
p_PrdDVI = plant.p(13 : end);

%% state equation and equilibrium constraints (rigid body)
% state equation (mass M and nonlinear term H)
M = M_func(plant, q);

C = C_func(plant, q, dq);
B = B_func();
W_N = W_N_func(plant, q);
W_T = W_T_func(plant, q);
H = -C + B' * tau_control + W_N' * [pN_1; pN_2; pN_3; pN_4] + W_T' * [pT_1; pT_2; pT_3; pT_4];

% equilibrium constraints
gap = gap_func(plant, q);

pT_1_lb = pT_1 - (-mu * pN_1);
pT_1_ub = (mu * pN_1) - pT_1;

pT_2_lb = pT_2 - (-mu * pN_2);
pT_2_ub = (mu * pN_2) - pT_2;

pT_3_lb = pT_3 - (-mu * pN_3);
pT_3_ub = (mu * pN_3) - pT_3;

pT_4_lb = pT_4 - (-mu * pN_4);
pT_4_ub = (mu * pN_4) - pT_4;

K_rigidBody = [gap(1);...
    pT_1_lb;...
    pT_1_ub;...
    gap(2);...
    pT_2_lb;...
    pT_2_ub;...
    gap(3);...
    pT_3_lb;...
    pT_3_ub;...
    gap(4);...
    pT_4_lb;...
    pT_4_ub];

l_rigidBody = zeros(4 * 3, 1);
u_rigidBody = Inf * ones(4 * 3, 1);

%% state equation and equilibrium constraints (periodic DVI)
% load parameter
alpha1 = plant.convergeSpeed1;
alpha2 = plant.convergeSpeed2;
mu_rlc = plant.radiusLimitCycle;
T = plant.periodLimitCycle;
beta = plant.dutyFactor;
b = plant.radiusSmoothVaries;
InitPhase = plant.InitPhase;
k_N = plant.k_N;
k_T = plant.k_T;

% connection matrix 
connectMatrix = cell(4, 4);
for i = 1 : 4
    for j = 1 : 4
        if i == j
            R_ij = zeros(2, 2); 
        else
            InitPhaseDiff_ij = InitPhase(i) - InitPhase(j);
            R_ij = [cos(InitPhaseDiff_ij), -sin(InitPhaseDiff_ij);...
                    sin(InitPhaseDiff_ij), cos(InitPhaseDiff_ij)];
        end
        connectMatrix{i, j} = R_ij;
    end
end

% state equation f and equilibrium constraints
f = sym('f', [plant.Dim.x - 2 * plant.qDim, 1]);
K_PrdDVI = sym('K', [4 * 2, 1]);
for i = 1 : 4
    % variable
    rho_i = x_PrdDVI(1 + (i - 1) * 4 : 2 + (i - 1) * 4);
    
    yN_i = p_PrdDVI(1 + (i - 1) * 2);
    yT_i = p_PrdDVI(2 + (i - 1) * 2); 
    
    % holf oscillator
    r_i = sqrt((rho_i(1))^2 + (rho_i(2))^2);
    gamma_i = pi / (beta * T * (exp(-b * rho_i(2)) + 1)) + ...
        pi / ((1 - beta) * T * (exp(b * rho_i(2)) + 1));
    f_HolfOsc_i = [alpha1 * (mu_rlc^2 - r_i^2) * rho_i(1) + gamma_i * rho_i(2);...
        alpha2 * (mu_rlc^2 - r_i^2) * rho_i(2) - gamma_i * rho_i(1)];
    sum_R_ij_rho_i = [0; 0]; % couple term
    for j = 1 : 4
        R_ij = connectMatrix{i, j};
        sum_R_ij_rho_i = sum_R_ij_rho_i + R_ij * rho_i;
    end
    
    % foot velocity trajectory 
    g_N_st_i = 0;
    g_N_sw_i = -k_N * sin(pi * rho_i(1)) * (f_HolfOsc_i(1) + sum_R_ij_rho_i(1));
    g_T_st_i = 0;
    g_T_sw_i = -k_T * cos(-pi/2 * rho_i(1)) * (f_HolfOsc_i(1) + sum_R_ij_rho_i(1));   
    f_fvt_i = [g_N_st_i * (1 - yN_i) + g_N_sw_i * yN_i;...
        g_T_st_i * (1 - yT_i) + g_T_sw_i * yT_i];
    
    % state equation
    f(1 + (i - 1) * 4 : 4 + (i - 1) * 4) = [f_HolfOsc_i + sum_R_ij_rho_i; f_fvt_i];   
    % equilibrium constraint
    K_PrdDVI(1 + (i - 1) * 2 : 2 + (i - 1) * 2) = [rho_i(2); rho_i(2)];    
end

l_PrdDVI = zeros(4 * 2, 1);
u_PrdDVI = ones(4 * 2, 1);

%% set dynamics
plant.M = M;
plant.H = H;
plant.f = f;

plant.K = [K_rigidBody; K_PrdDVI];

plant.l = [l_rigidBody; l_PrdDVI];
plant.u = [u_rigidBody; u_PrdDVI];

end

%% kinemaics and Jacobian
function kine = kinematics_1(plant, q, body, mode)
x = q(1);
z = q(2);
if strcmp(body, 'torso')
    l = plant.linkLength(1);
    d = plant.linkCenter(1);
elseif strcmp(body, 'thigh_1') || strcmp(body, 'thigh_2')
    l = plant.linkLength(2);
    d = plant.linkCenter(2);
end

switch body
    case 'torso'
        xita = q(3);        
    case 'thigh_1'
        xita = q(4);        
    case 'thigh_2'
        xita = q(6);               
end

switch mode
    case 'ee'
        kine = [x + l * sin(xita); z - l * cos(xita)];
    case 'com'
        kine = [x + d * sin(xita); z - d * cos(xita)];
end

end

function jac = jacobian_1(plant, q, body, mode)
if strcmp(body, 'torso')
    l = plant.linkLength(1);
    d = plant.linkCenter(1);
elseif strcmp(body, 'thigh_1') || strcmp(body, 'thigh_2')
    l = plant.linkLength(2);
    d = plant.linkCenter(2);
end

switch body
    case 'torso'
        xita = q(3);        
        switch mode
            case 'ee'
                jac = [1, 0, l * cos(xita), 0, 0, 0, 0, 0, 0, 0, 0;...
                       0, 1, l * sin(xita), 0, 0, 0, 0, 0, 0, 0, 0];
            case 'com'
                jac = [1, 0, d * cos(xita), 0, 0, 0, 0, 0, 0, 0, 0;...
                       0, 1, d * sin(xita), 0, 0, 0, 0, 0, 0, 0, 0];                
        end
    case 'thigh_1'
        xita = q(4); 
        switch mode
            case 'ee'
                jac = [1, 0, 0, l * cos(xita), 0, 0, 0, 0, 0, 0, 0;...
                       0, 1, 0, l * sin(xita), 0, 0, 0, 0, 0, 0, 0];
            case 'com'
                jac = [1, 0, 0, d * cos(xita), 0, 0, 0, 0, 0, 0, 0;...
                       0, 1, 0, d * sin(xita), 0, 0, 0, 0, 0, 0, 0];                
        end
    case 'thigh_2'
        xita = q(6);   
        switch mode
            case 'ee'
                jac = [1, 0, 0, 0, 0, l * cos(xita), 0, 0, 0, 0, 0;...
                       0, 1, 0, 0, 0, l * sin(xita), 0, 0, 0, 0, 0];
            case 'com'
                jac = [1, 0, 0, 0, 0, d * cos(xita), 0, 0, 0, 0, 0;...
                       0, 1, 0, 0, 0, d * sin(xita), 0, 0, 0, 0, 0];                
        end
end

end

function kine = kinematics_2(plant, q, body, mode)
if strcmp(body, 'calf_1') || strcmp(body, 'calf_2')
    l = plant.linkLength(3);
    d = plant.linkCenter(3);
elseif strcmp(body, 'thigh_3') || strcmp(body, 'thigh_4')
    l = plant.linkLength(2);
    d = plant.linkCenter(2);
end

switch body
    case 'calf_1'
        p = kinematics_1(plant, q, 'thigh_1', 'ee');
        xita = q(5);
    case 'calf_2'
        p = kinematics_1(plant, q, 'thigh_2', 'ee');
        xita = q(7);      
    case 'thigh_3'
        p = kinematics_1(plant, q, 'torso', 'ee');
        xita = q(8);
    case 'thigh_4'
        p = kinematics_1(plant, q, 'torso', 'ee');
        xita = q(10);
end

switch mode
    case 'ee'
        kine = p + [l * sin(xita); -l*cos(xita)];
    case 'com'
        kine = p + [d * sin(xita); -d*cos(xita)];
end

end

function jac = jacobian_2(plant, q, body, mode)
if strcmp(body, 'calf_1') || strcmp(body, 'calf_2')
    l = plant.linkLength(3);
    d = plant.linkCenter(3);
elseif strcmp(body, 'thigh_3') || strcmp(body, 'thigh_4')
    l = plant.linkLength(2);
    d = plant.linkCenter(2);
end

switch body
    case 'calf_1'
        jac = jacobian_1(plant, q, 'thigh_1', 'ee');
        xita = q(5);
        switch mode
            case 'ee'
                jac(1, 5) = jac(1, 5) + l * cos(xita);
                jac(2, 5) = jac(2, 5) + l * sin(xita);
            case 'com'
                jac(1, 5) = jac(1, 5) + d * cos(xita);
                jac(2, 5) = jac(2, 5) + d * sin(xita);                
        end
    case 'calf_2'
        jac = jacobian_1(plant, q, 'thigh_2', 'ee');
        xita = q(7);
        switch mode
            case 'ee'
                jac(1, 7) = jac(1, 7) + l * cos(xita);
                jac(2, 7) = jac(2, 7) + l * sin(xita);
            case 'com'
                jac(1, 7) = jac(1, 7) + d * cos(xita);
                jac(2, 7) = jac(2, 7) + d * sin(xita);                
        end       
    case 'thigh_3'
        jac = jacobian_1(plant, q, 'torso', 'ee');
        xita = q(8);
        switch mode
            case 'ee'
                jac(1, 8) = jac(1, 8) + l * cos(xita);
                jac(2, 8) = jac(2, 8) + l * sin(xita);
            case 'com'
                jac(1, 8) = jac(1, 8) + d * cos(xita);
                jac(2, 8) = jac(2, 8) + d * sin(xita);                
        end
    case 'thigh_4'
        jac = jacobian_1(plant, q, 'torso', 'ee');
        xita = q(10);
        switch mode
            case 'ee'
                jac(1, 10) = jac(1, 10) + l * cos(xita);
                jac(2, 10) = jac(2, 10) + l * sin(xita);
            case 'com'
                jac(1, 10) = jac(1, 10) + l * cos(xita);
                jac(2, 10) = jac(2, 10) + l * sin(xita);
        end
end

end

function kine = kinematics_3(plant, q, body, mode)
l = plant.linkLength(3);
d = plant.linkCenter(3);
switch body
    case 'calf_3'
        p = kinematics_2(plant, q, 'thigh_3', 'ee');
        xita = q(9);
    case 'calf_4'
        p = kinematics_2(plant, q, 'thigh_4', 'ee');
        xita = q(11);
end

switch mode
    case 'ee'
        kine = p + [l * sin(xita); -l * cos(xita)];
    case 'com'
        kine = p + [d * sin(xita); -d * cos(xita)];
end

end

function jac = jacobian_3(plant, q, body, mode)
l = plant.linkLength(3);
d = plant.linkCenter(3);
switch body
    case 'calf_3'
        jac = jacobian_2(plant, q, 'thigh_3', 'ee');
        xita = q(9);
        switch mode
            case 'ee'
                jac(1, 9) = jac(1, 9) + l * cos(xita);
                jac(2, 9) = jac(2, 9) + l * sin(xita);
            case 'com'
                jac(1, 9) = jac(1, 9) + d * cos(xita);
                jac(2, 9) = jac(2, 9) + d * sin(xita);
        end
    case 'calf_4'
        jac = jacobian_2(plant, q, 'thigh_4', 'ee');
        xita = q(11);
        switch mode
            case 'ee'
                jac(1, 11) = jac(1, 11) + l * cos(xita);
                jac(2, 11) = jac(2, 11) + l * sin(xita);
            case 'com'
                jac(1, 11) = jac(1, 11) + d * cos(xita);
                jac(2, 11) = jac(2, 11) + d * sin(xita);
        end
end

end

%% Lagrangian
function L = lagrangian(plant, q, dq)
m_torso = plant.mass(1);
m_thigh = plant.mass(2);
m_calf = plant.mass(3);
I_torso = plant.inertia(1);
I_thigh = plant.inertia(2);
I_calf = plant.inertia(3);
g = plant.g;

L = 0;

% torso
p_torso = kinematics_1(plant, q, 'torso', 'com');
J_torso = jacobian_1(plant, q, 'torso', 'com');
v_torso = J_torso * dq;

L = L + 0.5 * m_torso * v_torso' * v_torso;
L = L + 0.5 * I_torso * dq(3)^2;
L = L - m_torso * g * p_torso(2);

% thigh 1
p_thigh_1 = kinematics_1(plant, q, 'thigh_1', 'com');
J_thigh_1 = jacobian_1(plant, q, 'thigh_1', 'com');
v_thigh_1 = J_thigh_1 * dq;

L = L + 0.5 * m_thigh * v_thigh_1' * v_thigh_1;
L = L + 0.5 * I_thigh * dq(4)^2;
L = L - m_thigh * g * p_thigh_1(2);

% leg 1
p_calf_1 = kinematics_2(plant, q, 'calf_1', 'com');
J_calf_1 = jacobian_2(plant, q, 'calf_1', 'com');
v_calf_1 = J_calf_1 * dq;

L = L + 0.5 * m_calf * v_calf_1' * v_calf_1;
L = L + 0.5 * I_calf * dq(5)^2;
L = L - m_calf * g * p_calf_1(2);

% thigh 2
p_thigh_2 = kinematics_1(plant, q, 'thigh_2', 'com');
J_thigh_2 = jacobian_1(plant, q, 'thigh_2', 'com');
v_thigh_2 = J_thigh_2 * dq;

L = L + 0.5 * m_thigh * v_thigh_2' * v_thigh_2;
L = L + 0.5 * I_thigh * dq(6)^2;
L = L - m_thigh * g * p_thigh_2(2);

% leg 2
p_calf_2 = kinematics_2(plant, q, 'calf_2', 'com');
J_calf_2 = jacobian_2(plant, q, 'calf_2', 'com');
v_calf_2 = J_calf_2 * dq;

L = L + 0.5 * m_calf * v_calf_2' * v_calf_2;
L = L + 0.5 * I_calf * dq(7)^2;
L = L - m_calf * g * p_calf_2(2);

% thigh 3
p_thigh_3 = kinematics_2(plant, q, 'thigh_3', 'com');
J_thigh_3 = jacobian_2(plant, q, 'thigh_3', 'com');
v_thigh_3 = J_thigh_3 * dq;

L = L + 0.5 * m_thigh * v_thigh_3' * v_thigh_3;
L = L + 0.5 * I_thigh * dq(8)^2;
L = L - m_thigh * g * p_thigh_3(2);

% leg 3
p_calf_3 = kinematics_3(plant, q, 'calf_3', 'com');
J_calf_3 = jacobian_3(plant, q, 'calf_3', 'com');
v_calf_3 = J_calf_3 * dq;

L = L + 0.5 * m_calf * v_calf_3' * v_calf_3;
L = L + 0.5 * I_calf * dq(9)^2;
L = L - m_calf * g * p_calf_3(2);

% thigh 4
p_thigh_4 = kinematics_2(plant, q, 'thigh_4', 'com');
J_thigh_4 = jacobian_2(plant, q, 'thigh_4', 'com');
v_thigh_4 = J_thigh_4 * dq;

L = L + 0.5 * m_thigh * v_thigh_4' * v_thigh_4;
L = L + 0.5 * I_thigh * dq(10)^2;
L = L - m_thigh * g * p_thigh_4(2);

% leg 4
p_calf_4 = kinematics_3(plant, q, 'calf_4', 'com');
J_calf_4 = jacobian_3(plant, q, 'calf_4', 'com');
v_calf_4 = J_calf_4 * dq;

L = L + 0.5 * m_calf * v_calf_4' * v_calf_4;
L = L + 0.5 * I_calf * dq(11)^2;
L = L - m_calf * g * p_calf_4(2);

end

function dLdq = dLdq_func(plant, q, dq)

L = lagrangian(plant, q, dq);
dLdq = jacobian(L, q);

end

function dLddq = dLddq_func(plant, q, dq)

L = lagrangian(plant, q, dq);
dLddq = jacobian(L, dq);

end

%% dynamics matrix
function M = M_func(plant, q)
m_torso = plant.mass(1);
m_thigh = plant.mass(2);
m_calf = plant.mass(3);
I_torso = plant.inertia(1);
I_thigh = plant.inertia(2);
I_calf = plant.inertia(3);

M = diag([0, 0, I_torso, I_thigh, I_calf, I_thigh, I_calf, I_thigh, I_calf, I_thigh, I_calf]);
% torso
J_torso = jacobian_1(plant, q, 'torso', 'com');
M = M + m_torso * J_torso' * J_torso;

% thigh 1
J_thigh_1 = jacobian_1(plant, q, 'thigh_1', 'com');
M = M + m_thigh * J_thigh_1' * J_thigh_1;

% leg 1
J_calf_1 = jacobian_2(plant, q, 'calf_1', 'com');
M = M + m_calf * J_calf_1' * J_calf_1;

% thigh 2
J_thigh_2 = jacobian_1(plant, q, 'thigh_2', 'com');
M = M + m_thigh * J_thigh_2' * J_thigh_2;

% leg 2
J_calf_2 = jacobian_2(plant, q, 'calf_2', 'com');
M = M + m_calf * J_calf_2' * J_calf_2;

% thigh 3
J_thigh_3 = jacobian_2(plant, q, 'thigh_3', 'com');
M = M + m_thigh * J_thigh_3' * J_thigh_3;

% leg 3
J_calf_3 = jacobian_3(plant, q, 'calf_3', 'com');
M = M + m_calf * J_calf_3' * J_calf_3;

% thigh 4
J_thigh_4 = jacobian_2(plant, q, 'thigh_4', 'com');
M = M + m_thigh * J_thigh_4' * J_thigh_4;

% leg 4
J_calf_4 = jacobian_3(plant, q, 'calf_4', 'com');
M = M + m_calf * J_calf_4' * J_calf_4;

end

function C = C_func(plant, q, dq)

dLdq = dLdq_func(plant, q, dq);
dLddq = dLddq_func(plant, q, dq);

dLddq_dq = jacobian(dLddq, q);

C =  dLddq_dq*dq - dLdq';
end

function gap = gap_func(plant, q)
p_calf_1 = kinematics_2(plant, q, 'calf_1', 'ee');
p_calf_2 = kinematics_2(plant, q, 'calf_2', 'ee');
p_calf_3 = kinematics_3(plant, q, 'calf_3', 'ee');
p_calf_4 = kinematics_3(plant, q, 'calf_4', 'ee');

gap = [p_calf_1(2), p_calf_2(2), p_calf_3(2), p_calf_4(2)];
end

function B = B_func()
B = [0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0;...
     0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0;...
     0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0;...
     0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0;...
     0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0;...
     0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0;...
     0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0;...
     0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1];
end

function W_N = W_N_func(plant, q)
J_calf_1 = jacobian_2(plant, q, 'calf_1', 'ee');
J_calf_2 = jacobian_2(plant, q, 'calf_2', 'ee');
J_calf_3 = jacobian_3(plant, q, 'calf_3', 'ee');
J_calf_4 = jacobian_3(plant, q, 'calf_4', 'ee');

W_N = [J_calf_1(2, :);...
       J_calf_2(2, :);...
       J_calf_3(2, :);...
       J_calf_4(2, :)];
end

function W_T = W_T_func(plant, q)
J_calf_1 = jacobian_2(plant, q, 'calf_1', 'ee');
J_calf_2 = jacobian_2(plant, q, 'calf_2', 'ee');
J_calf_3 = jacobian_3(plant, q, 'calf_3', 'ee');
J_calf_4 = jacobian_3(plant, q, 'calf_4', 'ee');

W_T = [J_calf_1(1, :);...
       J_calf_2(1, :);...
       J_calf_3(1, :);...
       J_calf_4(1, :)];
end

