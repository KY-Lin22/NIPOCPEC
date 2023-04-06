function plotSimuResult(plant, simuTimeStep, InitState, tau, x, p)
%plotSimuResult
%   Detailed explanation goes here
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
plot(timeAxis, [InitState(1), x(1,:)], 'r',...
     timeAxis, [InitState(2), x(2,:)], 'g',...
     timeAxis, [InitState(3), x(3,:)], 'b',...
     timeAxis, [InitState(4), x(4,:)], 'k', 'LineWidth',1.2)
legend('x_b', 'y_b', '\theta_b', 'l_l')
legend('Orientation','horizontal')
xlabel('time [s]')
title('position')

subplot(3,1,2)
plot(timeAxis, [InitState(5), x(5,:)], 'r',...
     timeAxis, [InitState(6), x(6,:)], 'g',...
     timeAxis, [InitState(7), x(7,:)], 'b',...
     timeAxis, [InitState(8), x(8,:)], 'k', 'LineWidth',1.2)
xlabel('time [s]')
title('velocity')

subplot(3,1,3)
plot(timeAxis(2:end), tau(1,:),'b',...
    timeAxis(2:end), tau(2, :), 'k','LineWidth', 1.2)
legend('\tau_b', '\tau_l')
legend('Orientation','horizontal')
xlabel('time [s]')
title('control')

figure(112)
subplot(4,1,1)
plot(timeAxis(2:end), K(1, :), 'b','LineWidth', 1.2)
xlabel('time [s]')
ylabel('gap [m]')
title('normal impact contact')
subplot(4,1,2)
plot(timeAxis(2:end), p(1, :), 'k','LineWidth', 1.2)
xlabel('time [s]')
ylabel('impact [N]')

subplot(4,1,3)
plot(timeAxis(2:end), (p(2, :) - p(3, :)), 'b', 'LineWidth', 1.2)
xlabel('time [s]')
ylabel('foot vel [m/s]')
title('tangential friction contact')
subplot(4,1,4)
plot(timeAxis(2:end), tau(3, :), 'k', 'LineWidth', 1.2)
xlabel('time [s]')
ylabel('friction[N]')

end

