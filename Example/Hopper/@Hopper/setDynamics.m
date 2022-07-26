function setDynamics(plant)
%setDynamics
%   Detailed explanation goes here

q = plant.x(1 : 4);
dq = plant.x(5 : end);

pN = plant.p(1);
pT = plant.tau(3);

mb = plant.mass(1);
ml = plant.mass(2);
Ib = plant.inertia(1);
Il = plant.inertia(2);
mu = plant.mu;
g = plant.g;

%% set state equation
% set mass matrix
M = diag([mb + ml, mb + ml, Ib + Il, ml]);

% set nonlinear matrix
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
f = [dq;...
    inv(M)*H];
plant.f = f;

%% set equilibrium dynamics
Gap = q(2) - q(4) * cos(q(3));
pT_lb = pT - (-mu * pN);
pT_ub = (mu * pN) - pT;

K = [Gap;...
     pT_lb;...
     pT_ub];
plant.K = K;

l = [0;...
    0;...
    0];
plant.l = l;

u = [Inf;...
    Inf;...
    Inf];
plant.u = u;
end

