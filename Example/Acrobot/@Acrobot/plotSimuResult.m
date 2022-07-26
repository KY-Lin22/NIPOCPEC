function plotSimuResult(plant, simuTimeStep, InitState, tau, x, p)
%plotSimuResult
%   Detailed explanation goes here
Dim = plant.Dim;
nStages = size(x, 2);
timeAxis = 0 : simuTimeStep : nStages * simuTimeStep;

% compute value of K
K = zeros(Dim.p, nStages);
for n = 1 : nStages
    K(:, n) = plant.computeVIFunc(tau(:, n), x(:, n), p(:, n));
end

%% plot
figure(111)
subplot(2,3,1)
plot(timeAxis, [InitState(1), x(1,:)]*360/(2*pi), 'r',...
     timeAxis, [InitState(2), x(2,:)]*360/(2*pi), 'g','LineWidth',1.2)
legend('link 1', 'link 2')
xlabel('time(s)')
title('position(deg)')

subplot(2,3,2)
plot(timeAxis, [InitState(3), x(3,:)], 'r',...
     timeAxis, [InitState(4), x(4,:)], 'g', 'LineWidth',1.2)
legend('link 1', 'link 2')
xlabel('time(s)')
title('velocity')

subplot(2,3,3)
plot(timeAxis(2:end), tau(1,:), 'LineWidth', 1.2)
xlabel('time(s)')
title('control')

subplot(2,3,4)
plot(timeAxis(2:end), p(1, :), 'k',...
     timeAxis(2:end), K(1, :), 'b','LineWidth', 1.2)
legend('force', 'q2 - qmin') 
xlabel('time(s)')
title('equilibrium dynamics')

subplot(2,3,5)
plot(timeAxis(2:end), p(2, :), 'k',...
     timeAxis(2:end), K(2, :), 'b',...
     'LineWidth', 1.2)
legend('force', 'qmax - q2') 
xlabel('time(s)')
title('equilibrium dynamics')

end

