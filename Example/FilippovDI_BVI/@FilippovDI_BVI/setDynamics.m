function setDynamics(plant)
%setDynamics
%   Detailed explanation goes here

y = plant.p;
x = plant.x;
%% set state equation
f1 = 1; % switch function > 0
f2 = 3; % switch function < 0
f = f1 * (1 - y) + f2 * y;

plant.f = f;

%% set equilibrium dynamics
plant.K = x;

l = 0;
plant.l = l;

u = 1;
plant.u = u;
end

