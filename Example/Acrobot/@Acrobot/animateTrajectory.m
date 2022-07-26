function animateTrajectory(plant, simuTimeStep, InitState, tau, x, p)
%animateTrajectory
%   Detailed explanation goes here
linkLength = plant.linkLength;

% define time axis
nStages = size(x, 2);
timeAxis = 0 : simuTimeStep : nStages * simuTimeStep;

% get Acrobot configuration sequence based on given q sequence
Acrobot_X = zeros(8, nStages + 1);
Acrobot_Y = zeros(8, nStages + 1);
for n = 1 : nStages + 1
    if n == 1
        x_n = InitState;
    else
        x_n = x(:, n - 1);
    end
    [Acrobot_X(:, n), Acrobot_Y(:, n)] = getAcRobotConfiguration(plant, x_n(1:2, :));
end

% generate a trajectory-time sequence about contact force p for animation
p = [p(:, 1), p];
trajp1_X = cell(1,nStages + 1);
trajp1_Y = cell(1,nStages + 1);
trajp2_X = cell(1,nStages + 1);
trajp2_Y = cell(1,nStages + 1);
for n = 1 : nStages + 1
    trajp1_X{1, n} = [timeAxis(1:n), repmat(timeAxis(n), 1, nStages + 1 - n)];
    trajp1_Y{1, n} = [p(1, 1:n), repmat(p(1, n), 1, nStages + 1 - n)];
    trajp2_X{1, n} = [timeAxis(1:n), repmat(timeAxis(n), 1, nStages + 1 - n)];
    trajp2_Y{1, n} = [p(2, 1:n), repmat(p(2, n), 1, nStages + 1 - n)];
end

%%
% get figure size
figure(100)
pos = get(gcf, 'Position');
width = pos(3);
height = pos(4);
% define movie record
mov = zeros(height*3/2, width*3/2, 1, nStages + 1, 'uint8');
% pre allocate
subplot(3,2,[1, 2, 3, 4])
Acrobot = plot(Acrobot_X(1 : 4, 1), Acrobot_Y(1 : 4, 1), '.-k', 'MarkerSize', 10, 'LineWidth', 1.2);
hold on
maxLimit = plot(Acrobot_X(5 : 6, 1), Acrobot_Y(5 : 6, 1), '--b', 'MarkerSize', 10, 'LineWidth', 0.5);
hold on  
minLimit = plot(Acrobot_X(7 : 8, 1), Acrobot_Y(7 : 8, 1), '--b', 'MarkerSize', 10, 'LineWidth', 0.5);
axisLimit_X = [-1 * (linkLength(1) + linkLength(2));...
               1 * (linkLength(1) + linkLength(2))];
axisLimit_Y = [-1 * (linkLength(1) + linkLength(2));...
               1 * (linkLength(1) + linkLength(2))];
axis([axisLimit_X; axisLimit_Y]);
timePrint = title(sprintf('Time: %0.2f sec', timeAxis(1)));

subplot(3,2,[5, 6])
trajp1 = plot(trajp1_X{1, 1}, trajp1_Y{1, 1}, 'LineWidth', 1.2);
hold on
trajp2 = plot(trajp2_X{1, 1}, trajp2_Y{1, 1}, 'LineWidth', 1.2);
axis([0; timeAxis(end); min([p(1, :), p(2, :)]); max([p(1, :), p(2, :)])]);
legend('qmin', 'qmax')  
xlabel('time(s)')
ylabel('force(N)')
title('joint limit contact force')

%% animate trajectory
for n = 1 : nStages + 1
    % update XData and YData
    set(Acrobot, 'XData', Acrobot_X(1 : 4, n), 'YData', Acrobot_Y(1 : 4, n));
    set(maxLimit, 'XData', Acrobot_X(5 : 6, n), 'YData', Acrobot_Y(5 : 6, n));
    set(minLimit, 'XData', Acrobot_X(7 : 8, n), 'YData', Acrobot_Y(7 : 8, n));
    set(timePrint, 'String', sprintf('Time: %0.2f sec', timeAxis(n)));
    set(trajp1, 'XData', trajp1_X{1, n}, 'YData', trajp1_Y{1, n});
    set(trajp2, 'XData', trajp2_X{1, n}, 'YData', trajp2_Y{1, n});
    % get frame as an image
    f = getframe(gcf);
    % Create a colormap for the first frame. for the rest of the frames, use the same colormap
    if n == 1
        [mov(:,:,1,n), map] = rgb2ind(f.cdata, 256, 'nodither');
    else
        mov(:,:,1,n) = rgb2ind(f.cdata, map, 'nodither');
    end
end
% create an animated GIF
imwrite(mov, map, 'Acrobot.gif', 'DelayTime', 0, 'LoopCount', inf)

end

%%
function [AcRobot_X, AcRobot_Y] = getAcRobotConfiguration(plant, q)
L = plant.linkLength;
q2_max = plant.q2_max;
q2_min = plant.q2_min;
% link 1 coordinate
link_1_X = [0;...
            L(1) * sin(q(1))];
link_1_Y = [0;...
            -L(1) * cos(q(1))];
% link 2 coordinate
link_2_X = [L(1) * sin(q(1));...
            L(1) * sin(q(1)) + L(2) * sin(q(1) + q(2))];
link_2_Y = [-L(1) * cos(q(1));...
            -L(1) * cos(q(1)) - L(2) * cos(q(1) + q(2))];
%        
maxLimit_X = [L(1) * sin(q(1));...
              L(1) * sin(q(1)) + 0.5 * L(2) * sin(q(1) + q2_max)];   
maxLimit_Y = [-L(1) * cos(q(1));...
              -L(1) * cos(q(1)) - 0.5 * L(2) * cos(q(1) + q2_max)];          
          
minLimit_X = [L(1) * sin(q(1));...
              L(1) * sin(q(1)) + 0.5 * L(2) * sin(q(1) + q2_min)];   
minLimit_Y = [-L(1) * cos(q(1));...
              -L(1) * cos(q(1)) - 0.5 * L(2) * cos(q(1) + q2_min)]; 
% link coordinate
AcRobot_X = [link_1_X; link_2_X; maxLimit_X; minLimit_X];
AcRobot_Y = [link_1_Y; link_2_Y; maxLimit_Y; minLimit_Y];
end