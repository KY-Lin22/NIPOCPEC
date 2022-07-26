function setDynamics(plant)
%setDynamics
%   Detailed explanation goes here
%% set state equation
A = [1, -3; ...
    -8, 10];
B = [-3;...
    -1];
F = [4;...
    8];
f = A * plant.x + B * plant.p + F * plant.tau;

plant.f = f;

%% set equilibrium dynamics
C = [1, -3];
D = 5;
E = 3;
plant.K = C * plant.x + D * plant.p + E * plant.tau;

plant.l = -1;
plant.u = 1;

end

