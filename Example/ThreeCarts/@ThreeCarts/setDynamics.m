function setDynamics(plant)
%setDynamics
%   Detailed explanation goes here

q = plant.x(1 : 3, 1);
dq = plant.x(4 : end, 1);

%% set state equation
f = [dq;...
     (1/plant.mass(1)) * (plant.tau(1) - plant.viscousDamping * dq(1) - plant.p(1));...
     (1/plant.mass(2)) * (plant.tau(2) - plant.viscousDamping * dq(2) + plant.p(1) - plant.p(2));...
     (1/plant.mass(3)) * (plant.tau(3) - plant.viscousDamping * dq(3) + plant.p(2))];
plant.f = f;

%% set equilibrium dynamics
K = [q(2) - q(1) - 0.5 * plant.cartLength(2) - 0.5 * plant.cartLength(1);...
     q(3) - q(2) - 0.5 * plant.cartLength(3) - 0.5 * plant.cartLength(2)];
plant.K = K;

l = [0;...
     0];
plant.l = l;

u = [Inf;...
     Inf];
plant.u = u;

end

