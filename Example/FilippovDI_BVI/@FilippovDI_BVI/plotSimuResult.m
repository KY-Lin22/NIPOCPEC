function plotSimuResult(plant, simuTimeStep, InitState, tau, x, p)
%plotSimuResult
%   Detailed explanation goes here
%%
Dim = plant.Dim;
nStages = size(x, 2);
timeAxis = 0 : simuTimeStep : nStages * simuTimeStep;

% y and z
y = p(1, :);
z = x(1, :);

% f
f = zeros(Dim.x, nStages);
for n = 1 : nStages
    f(:, n) = plant.computeStateEquation_Function(tau(:, n), x(:, n), p(:, n));
end

%% plot
figure(112)
subplot(3,2,1)
plot(timeAxis, [InitState, x], 'r', 'LineWidth',1.2)
legend('x')
xlabel('time(s)')
title('system state')

subplot(3,2,2)
plot(timeAxis(2:end), f, 'b', 'LineWidth', 1.2)
legend('f') 
xlabel('time(s)')
title('state equation')

subplot(3,2,3)
plot(timeAxis(2:end), y, 'k', 'LineWidth', 1.2)
legend('y') 
xlabel('time(s)')
title('smoothing function')

subplot(3,2,4)
plot(timeAxis(2:end), z, 'g', 'LineWidth', 1.2)
legend('z := x') 
xlabel('time(s)')
title('switch function')

subplot(3,2,5)
plot(timeAxis(2:end), z, 'k',...
    timeAxis(2:end), y, 'b', 'LineWidth', 1.2)
legend('z(K)', 'y(p)') 
xlabel('time(s)')
title('checking BVI')

end

