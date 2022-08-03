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
subplot(3,2,1)
plot(timeAxis, [InitState(1), x(1, :)], 'r',...
     timeAxis, [InitState(2), x(2, :)], 'g',...
     timeAxis, [InitState(3), x(3, :)], 'b', 'LineWidth', 1.2)
legend('cart 1', 'cart 2', 'cart 3')
xlabel('time(s)')
ylabel('position (m)')
title('cart position')

subplot(3,2,3)
plot(timeAxis, [InitState(4), x(4, :)], 'r',...
     timeAxis, [InitState(5), x(5, :)], 'g',...
     timeAxis, [InitState(6), x(6, :)], 'b', 'LineWidth', 1.2)
legend('cart 1', 'cart 2', 'cart 3')
xlabel('time(s)')
ylabel('rate(m/s)')
title('cart velocity')

subplot(3,2,5)
plot(timeAxis(2:end), tau(1,:), 'r',...
     timeAxis(2:end), tau(2,:), 'g',...
     timeAxis(2:end), tau(3,:), 'b', 'LineWidth', 1.2)
legend('cart 1', 'cart 2', 'cart 3')
xlabel('time(s)')
ylabel('torque(Nm)')
title('actuated force')

subplot(3,2,2)
plot(timeAxis(2:end), K(1, :),...
     timeAxis(2:end), K(2, :), 'LineWidth', 1.2)
legend('cart 2 -> 1', 'cart 3 -> 2') 
xlabel('time(s)')
ylabel('distance(m)')
title('gap')

subplot(3,2,4)
plot(timeAxis(2:end), p(1,:),...
     timeAxis(2:end), p(2,:),'LineWidth', 1.2)
legend('cart 2 -> 1', 'cart 3 -> 2') 
xlabel('time(s)')
ylabel('force(N)')
title('normal contact force')

end

