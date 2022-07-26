function setDynamics(plant)
%setDynamics
%   Detailed explanation goes here

q = plant.x(1 : 2, 1);
qC = q(1);
qP = q(2);
dq = plant.x(3 : end, 1);
dqC = dq(1);
dqP = dq(2);
mass = plant.mass;
linkLength = plant.linkLength;
g = plant.g;

%% set state equation
% set mass matrix
M = [mass(1) + mass(2),                    mass(2) * linkLength(1) * cos(qP);...
     mass(2) * linkLength(1) * cos(qP),  mass(2) * linkLength(1)^2];

% set nonlinear matrix including control torque        
C = [0,   -mass(2) * linkLength(1) * dqP * sin(qP);...
     0,   0];

G = [0;...
     -mass(2) * g * linkLength(1) * sin(qP)];

Bu = [plant.tau(1);...
      0]; 
P = [plant.p(1);...
     0]; % friction bewteen cart and ground

H = G + Bu + P - C * [dqC; dqP];   

f = [dq;...
     inv(M)*H];
plant.f = f;

%% set equilibrium dynamics
plant.K = dqC;

plant.l = -2;
plant.u = 2;

end

