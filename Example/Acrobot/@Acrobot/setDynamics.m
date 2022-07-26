function setDynamics(plant)
%setDynamics
%   Detailed explanation goes here

q = plant.x(1 : 2, 1);
dq = plant.x(3 : end, 1);

m = plant.mass;
L = plant.linkLength;
Lc = plant.linkCenter;
I = plant.inertia;
jf = plant.jointFriction;
g = plant.g;

%% set state equation
% set mass matrix
M = [I(1) + I(2) + m(2)*L(1)^2 + 2*m(2)*L(1)*Lc(2)*cos(q(2)), I(2) + m(2)*L(1)*Lc(2)*cos(q(2));...
     I(2) + m(2)*L(1)*Lc(2)*cos(q(2)),                        I(2)];

% set nonlinear matrix including control torque 
C = [-2*m(2)*L(1)*Lc(2)*sin(q(2))*dq(2),    -m(2)*L(1)*Lc(2)*sin(q(2))*dq(2);...
     m(2)*L(1)*Lc(2)*sin(q(2))*dq(1),        0];

G = [-m(1)*g*Lc(1)*sin(q(1)) - m(2)*g*(L(1)*sin(q(1)) + Lc(2)*sin(q(1)+q(2)));...
     -m(2)*g*Lc(2)*sin(q(1)+q(2))];
 
Bu = [0;...
      plant.tau(1)];
  
P = [0;...
     plant.p(1) - plant.p(2)];  

H = G + Bu + P - C * [dq(1); dq(2)] - [jf(1) * dq(1); jf(2) * dq(2)];   

f = [dq;...
     inv(M)*H];

plant.f = f;

%% set equilibrium dynamics
K = [q(2) - plant.q2_min;...
     plant.q2_max - q(2)];
plant.K = K;

l = [0;...
     0];
plant.l = l;

u = [Inf;...
     Inf];
plant.u = u;

end

