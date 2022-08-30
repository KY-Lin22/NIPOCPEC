function setDynamics(plant)
%setDynamics
%   Detailed explanation goes here

y = plant.tau(1);
control = plant.tau(2);

%% set state equation
f1 = 3;
f2 = 1;
f = f1 * (1 - y) + f2 * y + control;

plant.f = f;

%% set equilibrium dynamics
K0 = y;
K1 = 1 - y;
plant.K = [K0;...
    K1];

l = [0;...
    0];
plant.l = l;

u = [Inf;...
    Inf];
plant.u = u;
end
