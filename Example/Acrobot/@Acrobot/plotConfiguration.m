function plotConfiguration(plant, InitState, RefState)
%plotConfiguration
%   Detailed explanation goes here

linkLength = plant.linkLength;

[initAcRobot_X, initAcRobot_Y] = getAcRobotConfiguration(plant, InitState(1 : 2, 1));
[refAcRobot_X, refAcRobot_Y] = getAcRobotConfiguration(plant, RefState(1 : 2, 1));

figure(2000);
subplot(1,2,1)
% plot InitPosition
plot(initAcRobot_X(1 : 4, 1), initAcRobot_Y(1 : 4, 1), '.-k',  'MarkerSize', 10, 'LineWidth', 1.2);% acrobot
hold on
plot(initAcRobot_X(5 : 6, 1), initAcRobot_Y(5 : 6, 1), '--b', 'MarkerSize', 10, 'LineWidth', 0.5);% maxLimit
hold on
plot(initAcRobot_X(7 : 8, 1), initAcRobot_Y(7 : 8, 1), '--b', 'MarkerSize', 10, 'LineWidth', 0.5);% minLimit
title('InitPosition')
axis equal
axisLimit_X = [-1.2 * (linkLength(1) + linkLength(2));...
              1.2 * (linkLength(1) + linkLength(2))];
axisLimit_Y = [-1.2 * (linkLength(1) + linkLength(2));...
              1.2 * (linkLength(1) + linkLength(2))];
axis([axisLimit_X; axisLimit_Y]);
% plot RefPosition
subplot(1,2,2)
plot(refAcRobot_X(1 : 4, 1), refAcRobot_Y(1 : 4, 1), '.-k',  'MarkerSize', 10, 'LineWidth', 1.2);% acrobot
hold on
plot(refAcRobot_X(5 : 6, 1), refAcRobot_Y(5 : 6, 1), '--b', 'MarkerSize', 10, 'LineWidth', 0.5);% maxLimit
hold on
plot(refAcRobot_X(7 : 8, 1), refAcRobot_Y(7 : 8, 1), '--b', 'MarkerSize', 10, 'LineWidth', 0.5);% minLimit
title('Ref Position')
axis equal
axisLimit_X = [-1.2 * (linkLength(1) + linkLength(2));...
              1.2 * (linkLength(1) + linkLength(2))];
axisLimit_Y = [-1.2 * (linkLength(1) + linkLength(2));...
              1.2 * (linkLength(1) + linkLength(2))];
axis([axisLimit_X; axisLimit_Y]);

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
