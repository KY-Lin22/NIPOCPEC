function setDynamics(plant)
%setDynamics
%   Detailed explanation goes here

q = plant.x(1 : 4);
dq = plant.x(5 : 8);
rho = plant.x(9 : 10);

pN = plant.p(1);
pT = plant.tau(3);

yT = plant.p(4);
yN = plant.p(5);

mb = plant.mass(1);
ml = plant.mass(2);
Ib = plant.inertia(1);
Il = plant.inertia(2);
mu_fc = plant.frictionCoeff;
g = plant.g;

a1 = plant.convergeSpeed1;
a2 = plant.convergeSpeed2;
mu_rlc = plant.radiusLimitCycle;
T = plant.periodLimitCycle;
beta = plant.dutyFactor;
b = plant.radiusSmoothVaries;

%% set state equation
% rigid body dynamics
M = diag([mb + ml, mb + ml, Ib + Il, ml]);
C = [0;...
    (mb + ml)*g;...
    0;...
    0];% here C := C*dq - G
B = [0, -sin(q(3));...
    0, cos(q(3));...
    1, 0;...
    0, 1];
W_N = [0, 1, q(4)*sin(q(3)), -cos(q(3))];
W_T = [1, 0, q(4)*cos(q(3)), sin(q(3))];

H = -C + B * [plant.tau(1); plant.tau(2)] + W_N' * pN + W_T' * pT;

f_rbd = inv(M)*H;

% holf oscillator
r = sqrt((rho(1))^2 + (rho(2))^2);
gamma = pi / (beta * T * (exp(-b * rho(2)) + 1)) + ...
    pi / ((1 - beta) * T * (exp(b * rho(2)) + 1));
f_HolfOsc = [a1 * (mu_rlc^2 - r^2) * rho(1) + gamma * rho(2);...
    a2 * (mu_rlc^2 - r^2) * rho(2) - gamma * rho(1)];

% foot velocity trajectory
k_N = 0.5;
k_T = 0.5;
g_N_st = 0;
g_N_sw = k_N * sin(pi * rho(1));
g_T_st = 0;
g_T_sw = k_T * cos((pi/2) * rho(1));

f_fvt = [g_N_st * (1 - yN) + g_N_sw * yN;...
    g_T_st * (1 - yT) + g_T_sw * yT];

%
f = [dq;...
    f_rbd;...
    f_HolfOsc;...
    f_fvt];
plant.f = f;


%% set equilibrium dynamics
Gap = q(2) - q(4) * cos(q(3));
pT_lb = pT - (-mu_fc * pN);
pT_ub = (mu_fc * pN) - pT;

K = [Gap;...
     pT_lb;...
     pT_ub;...
     rho(2);...
     rho(2)];
plant.K = K;

l = [0;...
    0;...
    0;...
    0;...
    0];
plant.l = l;

u = [Inf;...
    Inf;...
    Inf;...
    1;...
    1];
plant.u = u;


end

