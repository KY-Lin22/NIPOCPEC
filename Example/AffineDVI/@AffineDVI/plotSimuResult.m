function plotSimuResult(plant, simuTimeStep, InitState, tau, x, p)
%plotSimuResult
%   Detailed explanation goes here

%% 
Dim = plant.Dim;
nStages = size(x, 2);
timeAxis = 0 : simuTimeStep : nStages * simuTimeStep;

% compute value of K
K = zeros(Dim.p, nStages);
for n = 1 : nStages
    K(:, n) = plant.computeVI_Function(tau(:, n), x(:, n), p(:, n));
end

%% plot
figure(111)
subplot(3,1,1)
plot(timeAxis, [InitState(1), x(1, :)], 'r',...
     timeAxis, [InitState(2), x(2, :)], 'g', 'LineWidth',1.2)
legend('x1', 'x2')
xlabel('time(s)')
title('system state')

subplot(3,1,2)
plot(timeAxis(2:end), tau(1,:), 'LineWidth', 1.2)
xlabel('time(s)')
title('control')

subplot(3,1,3)
plot(timeAxis(2:end), p(1, :), 'k',...
     timeAxis(2:end), K(1, :), 'b', 'LineWidth', 1.2)
legend('p', 'K') 
xlabel('time(s)')
title('equilibrium dynamics')
end

