function setDynamics(plant)
%setDynamics
%   Detailed explanation goes here

q = plant.x(1 : 7, :);
dq = plant.x(8 : end, :);
tau_control = plant.tau(1 : 5, :);% control
pN_1 = plant.p(1);% normal contact force in leg 1
pT_1 = plant.tau(6); % friction force in leg 1

pN_2 = plant.p(4); % normal contact force in leg 2
pT_2 = plant.tau(7); % friction force in leg 2

mu = plant.mu;

%% set state equation
M = M_func(plant, q);

C = C_func(plant, q, dq);

B = B_func();

W_N = W_N_func(plant, q);

W_T = W_T_func(plant, q);

H = -C + B' * tau_control + W_N' * [pN_1; pN_2] + W_T' * [pT_1; pT_2];

plant.M = M;
plant.H = H;

%% set equilibrium dynamics
gap = gap_func(plant, q);

pT_1_lb = pT_1 - (-mu * pN_1);
pT_1_ub = (mu * pN_1) - pT_1;

pT_2_lb = pT_2 - (-mu * pN_2);
pT_2_ub = (mu * pN_2) - pT_2;

K = [gap(1);...
    pT_1_lb;...
    pT_1_ub;...
    gap(2);...
    pT_2_lb;...
    pT_2_ub];
plant.K = K;

plant.l = [0; 0; 0; 0; 0; 0];
plant.u = [Inf; Inf; Inf; Inf; Inf; Inf];

end

%% kinematics and Jacobian
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
        switch mode
            case 'ee'
                kine = [x - l * sin(xita); z + l * cos(xita)];
            case 'com'
                kine = [x - d * sin(xita); z + d * cos(xita)];
        end       
    case 'thigh_1'
        xita = q(4);
        switch mode
            case 'ee'
                kine = [x + l * sin(xita); z - l * cos(xita)];
            case 'com'
                kine = [x + d * sin(xita); z - d * cos(xita)];
        end
    case 'thigh_2'
        xita = q(6);   
        switch mode
            case 'ee'
                kine = [x + l * sin(xita); z - l * cos(xita)];
            case 'com'
                kine = [x + d * sin(xita); z - d * cos(xita)];
        end        
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
                jac = [1, 0, -l * cos(xita), 0, 0, 0, 0;...
                       0, 1, -l * sin(xita), 0, 0, 0, 0];            
            case 'com'
                jac = [1, 0, -d * cos(xita), 0, 0, 0, 0;...
                       0, 1, -d * sin(xita), 0, 0, 0, 0];
        end 
    case 'thigh_1'
        xita = q(4);
        switch mode
            case 'ee'
                jac = [1, 0, 0, l * cos(xita), 0, 0, 0;...
                       0, 1, 0, l * sin(xita), 0, 0, 0];
            case 'com'
                jac = [1, 0, 0, d * cos(xita), 0, 0, 0;...
                       0, 1, 0, d * sin(xita), 0, 0, 0];                
        end        
    case 'thigh_2'
        xita = q(6);   
        switch mode
            case 'ee'
                jac = [1, 0, 0, 0, 0, l * cos(xita), 0;...
                       0, 1, 0, 0, 0, l * sin(xita), 0];
            case 'com'
                jac = [1, 0, 0, 0, 0, d * cos(xita), 0;...
                       0, 1, 0, 0, 0, d * sin(xita), 0];
        end         
end

end


function kine = kinematics_2(plant, q, body, mode)
lb = plant.linkLength(3);
db = plant.linkCenter(3);
switch body
    case 'calf_1'
        p = kinematics_1(plant, q, 'thigh_1', 'ee');
        xita_b = q(5);
    case 'calf_2'
        p = kinematics_1(plant, q, 'thigh_2', 'ee');
        xita_b = q(7);
end
switch mode
    case 'ee'
        kine = p + [lb * sin(xita_b); -lb * cos(xita_b)];
    case 'com'
        kine = p + [db * sin(xita_b); -db * cos(xita_b)];
end

end


function jac = jacobian_2(plant, q, body, mode)
lb = plant.linkLength(3);
db = plant.linkCenter(3);
switch body
    case 'calf_1'
        jac = jacobian_1(plant, q, 'thigh_1', 'ee');
        xita_b = q(5);
        switch mode
            case 'ee'
                jac(1, 5) = jac(1, 5) + lb * cos(xita_b);
                jac(2, 5) = jac(2, 5) + lb * sin(xita_b);
            case 'com'
                jac(1, 5) = jac(1, 5) + db * cos(xita_b);
                jac(2, 5) = jac(2, 5) + db * sin(xita_b);
        end
    case 'calf_2'
        jac = jacobian_1(plant, q, 'thigh_2', 'ee');
        xita_b = q(7);
        switch mode
            case 'ee'
                jac(1, 7) = jac(1, 7) + lb * cos(xita_b);
                jac(2, 7) = jac(2, 7) + lb * sin(xita_b);               
            case 'com'
                jac(1, 7) = jac(1, 7) + db * cos(xita_b);
                jac(2, 7) = jac(2, 7) + db * sin(xita_b);
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

M = diag([0, 0, I_torso, I_thigh, I_calf, I_thigh, I_calf]);
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

gap = [p_calf_1(2), p_calf_2(2)];
end

function B = B_func()
B = [0, 0 ,1, 0, 0, 0, 0;...
    0, 0, 0, 1, 0, 0, 0;...
    0, 0, 0, 0, 1, 0, 0;...
    0, 0, 0, 0, 0, 1, 0;...
    0, 0, 0, 0, 0, 0, 1];
end

function W_N = W_N_func(plant, q)
J_calf_1 = jacobian_2(plant, q, 'calf_1', 'ee');
J_calf_2 = jacobian_2(plant, q, 'calf_2', 'ee');

W_N = [J_calf_1(2, :);...
    J_calf_2(2, :)];
end

function W_T = W_T_func(plant, q)
J_calf_1 = jacobian_2(plant, q, 'calf_1', 'ee');
J_calf_2 = jacobian_2(plant, q, 'calf_2', 'ee');

W_T = [J_calf_1(1, :);...
    J_calf_2(1, :)];
end