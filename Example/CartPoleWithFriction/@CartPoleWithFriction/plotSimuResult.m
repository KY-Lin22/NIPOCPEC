function plotSimuResult(plant, simuTimeStep, InitState, tau, x, p)
%plotSimuResult
%   Detailed explanation goes here
nStages = size(x, 2);
timeAxis = 0 : simuTimeStep : nStages * simuTimeStep;
Dim = plant.Dim;

K = zeros(Dim.p, nStages);
for n = 1 : nStages
    K(:, n) = plant.computeVI_Function(tau(:, n), x(:, n), p(:, n));
end

%% plot
figure(111)
subplot(4,1,1)
plot(timeAxis, [InitState(1), x(1, :)], 'r',...
     timeAxis, [InitState(2), x(2, :)], 'g','LineWidth',1.2)
legend('cart', 'pole')
xlabel('time [s]')
title('position')

subplot(4,1,2)
plot(timeAxis, [InitState(3), x(3, :)], 'r',...
     timeAxis, [InitState(4), x(4, :)], 'g', 'LineWidth',1.2)
xlabel('time [s]')
title('velocity')

subplot(4,1,3)
plot(timeAxis(2:end), tau(1,:), 'LineWidth', 1.2)
xlabel('time [s]')
title('control')

subplot(4,1,4)
plot(timeAxis(2:end), p(1, :), 'k',...
     timeAxis(2:end), K(1, :), 'b', 'LineWidth', 1.2)
legend('friction', 'cart vel') 
xlabel('time [s]')
title('equilibrium dynamics')

end

