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
subplot(3,3,1)
plot(timeAxis, [InitState(1), x(1,:)], 'r',...
     timeAxis, [InitState(2), x(2,:)], 'g',...
     timeAxis, [InitState(3), x(3,:)], 'b',...
     timeAxis, [InitState(4), x(4,:)], 'c',...
     timeAxis, [InitState(5), x(5,:)], 'm',...
     timeAxis, [InitState(6), x(6,:)], 'y',...
     timeAxis, [InitState(7), x(7,:)], 'k',...
     'LineWidth',1.2)
legend('x', 'z', 'torso', 'thigh1', 'calf1', 'thigh2', 'calf2')
xlabel('time(s)')
title('position')

subplot(3,3,2)
plot(timeAxis, [InitState(8), x(8,:)], 'r',...
     timeAxis, [InitState(9), x(9,:)], 'g',...
     timeAxis, [InitState(10), x(10,:)], 'b',...
     timeAxis, [InitState(11), x(11,:)], 'c',...
     timeAxis, [InitState(12), x(12,:)], 'm',...
     timeAxis, [InitState(13), x(13,:)], 'y',...
     timeAxis, [InitState(14), x(14,:)], 'k',...
     'LineWidth',1.2)
legend('x', 'z', 'torso', 'thigh1', 'calf1', 'thigh2', 'calf2')
xlabel('time(s)')
title('velocity')

subplot(3,3,3)
plot(timeAxis(2:end), tau(1,:),'b',...
    timeAxis(2:end), tau(2, :), 'c',...
    timeAxis(2:end), tau(3, :), 'm',...
    timeAxis(2:end), tau(4, :), 'y',...
    timeAxis(2:end), tau(5, :), 'k',...
    'LineWidth', 1.2)
legend('torso', 'thigh1', 'calf1', 'thigh2', 'calf2')
xlabel('time(s)')
title('control')

subplot(3,3,4)
plot(timeAxis(2:end), p(1, :), 'k',...
     timeAxis(2:end), K(1, :), 'b','LineWidth', 1.2)
legend('force', 'gap') 
xlabel('time(s)')
title('normal contact force (leg 1)')

subplot(3,3,5)
plot(timeAxis(2:end), p(4, :), 'k',...
     timeAxis(2:end), K(4, :), 'b','LineWidth', 1.2)
legend('force', 'gap') 
xlabel('time(s)')
title('normal contact force (leg 2)')

subplot(3,3,7)
plot(timeAxis(2:end), tau(6, :), 'k',...
     timeAxis(2:end), (p(2, :) - p(3, :)), 'b',...
     'LineWidth', 1.2)
legend('force', 'velocity') 
xlabel('time(s)')
title('friction (leg 1)')

subplot(3,3,8)
plot(timeAxis(2:end), tau(7, :), 'k',...
     timeAxis(2:end), (p(5, :) - p(6, :)), 'b',...
     'LineWidth', 1.2)
legend('force', 'velocity') 
xlabel('time(s)')
title('friction (leg 2)')

end

