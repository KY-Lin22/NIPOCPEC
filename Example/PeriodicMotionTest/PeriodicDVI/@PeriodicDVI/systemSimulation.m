function [state, phase, radius] = systemSimulation(plant, InitState, nStages, timeStep)
%simulatePeriodicDVI
%   Detailed explanation goes here
%% initialization
InitPhase = plant.InitPhase;
num_Osci = length(InitPhase);

state = zeros(4 * num_Osci, nStages);
phase = zeros(1 * num_Osci, nStages);
radius = zeros(1 * num_Osci, nStages);

%% simulation
for n = 1 : nStages
    % previous state (x: previous, x_n: current)
    if n == 1
        x = InitState;
    else
        x = state(:, n - 1);
    end
    % determine the value of auxiliary variable y by switch function rho(2)
    p = zeros(plant.Dim.p, 1);
    for i = 1 : num_Osci
        if x(2 + (i - 1) * 4) > 0
            p(1 + (i - 1) * 2 : 2 + (i - 1) * 2) = [0; 0];
        else
            p(1 + (i - 1) * 2 : 2 + (i - 1) * 2) = [1; 1];
        end
    end  
    %
    tau = zeros(plant.Dim.tau, 1);
    % Periodic DVI dynamics
    f = plant.computeStateEquation_Function(tau, x, p);
    % update current state (forward Euler method)
    x_n = x + f * timeStep;
    state(:, n) = x_n;    
    % update phase and radius
    for i = 1 : num_Osci
        rho_i_n = x_n(1 + (i - 1) * 4 : 2 + (i - 1) * 4);
        phase(i, n) = atan2(rho_i_n(2), rho_i_n(1));
        radius(i, n) = norm(rho_i_n, 2);
    end
end

end

