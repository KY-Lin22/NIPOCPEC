function plotSimuResult(plant, simuTimeStep, InitState, tau, x, p)
%plotSimuResult
%   Detailed explanation goes here
%%
Dim = plant.Dim;
nStages = size(x, 2);
timeAxis = 0 : simuTimeStep : nStages * simuTimeStep;

% y and z
y = tau(1, :);
z = p(2, :) - p(1, :);

% f
f = zeros(Dim.x, nStages);
for n = 1 : nStages
    f(:, n) = plant.computeStateEquation_Function(tau(:, n), x(:, n), p(:, n));
end

%% plot
figure(111)
subplot(4,1,1)
plot(timeAxis, [InitState, x], 'r', 'LineWidth',1.2)
legend('x')
xlabel('time(s)')
title('system state')

subplot(4,1,2)
plot(timeAxis(2:end), f, 'b', 'LineWidth', 1.2)
legend('f') 
xlabel('time(s)')
title('state equation')

subplot(4,1,3)
plot(timeAxis(2:end), y, 'k', 'LineWidth', 1.2)
legend('y') 
xlabel('time(s)')
title('smoothing function')

subplot(4,1,4)
plot(timeAxis(2:end), z, 'g', 'LineWidth', 1.2)
legend('z') 
xlabel('time(s)')
title('switch function')

figure(112)
subplot(2,1,1)
plot(timeAxis(2:end), p(1, :), 'k',...
    timeAxis(2:end), y, 'b', 'LineWidth', 1.2)
legend('lambda0', 'y') 
xlabel('time(s)')
title('NCP (lambda_0 and y)')

subplot(2,1,2)
plot(timeAxis(2:end), p(2, :), 'k',...
    timeAxis(2:end), ones(1, nStages) - y, 'b', 'LineWidth', 1.2)
legend('lambda1', '1-y') 
xlabel('time(s)')
title('NCP (lambda_1 and 1 - y)')
end

