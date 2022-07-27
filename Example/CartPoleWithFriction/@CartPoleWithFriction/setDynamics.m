function setDynamics(plant)
%setDynamics
%   Detailed explanation goes here

q = plant.x(1 : 2, 1);
dq = plant.x(3 : end, 1);
mass = plant.mass;
linkLength = plant.linkLength;
g = plant.g;

%% set state equation
% set mass matrix
M = [mass(1) + mass(2),                    mass(2) * linkLength * cos(q(2));...
     mass(2) * linkLength * cos(q(2)),  mass(2) * linkLength^2];

% set nonlinear matrix including control torque        
C = [0,   -mass(2) * linkLength * dq(2) * sin(q(2));...
     0,   0];

G = [0;...
     -mass(2) * g * linkLength * sin(q(2))];

Bu = [plant.tau(1);...
      0]; 
P = [plant.p(1);...
     0]; % friction bewteen cart and ground

H = G + Bu + P - C * [dq(1); dq(2)];   

f = [dq;...
     inv(M)*H];
plant.f = f;

%% set equilibrium dynamics
plant.K = dq(1);

plant.l = -2;
plant.u = 2;

end

