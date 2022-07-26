function plotConfiguration(plant, InitState, RefState)
%plotConfiguration
%   Detailed explanation goes here

DynVarLimit = plant.DynVarLimit;
%
baseline_X = [DynVarLimit.x_Min(1, 1) - 0.2; DynVarLimit.x_Max(1, 1) + 0.1];  
baseline_Y = [0;0];

[initQuadruped_X, initQuadruped_Y] = getQuadrupedConfiguration(plant, InitState(1:11, :));

refStageNum = size(RefState, 2);
refQuadruped_X = zeros(14, refStageNum);
refQuadruped_Y = zeros(14, refStageNum);
for n = 1 : refStageNum
    [refQuadruped_X(:, n), refQuadruped_Y(:, n)] = getQuadrupedConfiguration(plant, RefState(1:11, n));
end

figure(2000)
% plot init position
subplot(2, 2, [1, 2])
plot(baseline_X, baseline_Y, '-k', 'MarkerSize', 1, 'LineWidth', 2);
hold on
plot(initQuadruped_X(1:2, 1), initQuadruped_Y(1:2, 1), '.-k', 'MarkerSize', 15, 'LineWidth', 6)% torso
hold on
plot(initQuadruped_X(3:5, 1), initQuadruped_Y(3:5, 1), '.-', 'Color', [0.8500 0.3250 0.0980], 'MarkerSize', 15, 'LineWidth', 1)% leg 1
hold on
plot(initQuadruped_X(6:8, 1), initQuadruped_Y(6:8, 1), '.-', 'Color', [0 0.4470 0.7410], 'MarkerSize', 15, 'LineWidth', 1)% leg 2
hold on
plot(initQuadruped_X(9:11, 1), initQuadruped_Y(9:11, 1), '.-', 'Color', [0.8500 0.3250 0.0980], 'MarkerSize', 15, 'LineWidth', 1)% leg 3
hold on
plot(initQuadruped_X(12:14, 1), initQuadruped_Y(12:14, 1), '.-', 'Color', [0 0.4470 0.7410], 'MarkerSize', 15, 'LineWidth', 1)% leg 4

title('Init Position')
axis equal
axisLimit_X = baseline_X;
axisLimit_Y = [-0.1; 1.1];
axis([axisLimit_X; axisLimit_Y]);

% plot RefPosition
subplot(2,2, [3, 4])
plot(baseline_X, baseline_Y, '-k', 'MarkerSize', 1, 'LineWidth', 2);
for n = 1 : refStageNum
    hold on
    plot(refQuadruped_X(1:2, n), refQuadruped_Y(1:2, n), '.-k', 'MarkerSize', 15, 'LineWidth', 6)% torso
    hold on
    plot(refQuadruped_X(3:5, n), refQuadruped_Y(3:5, n), '.-', 'Color', [0.8500 0.3250 0.0980], 'MarkerSize', 15, 'LineWidth', 2)% leg 1
    hold on
    plot(refQuadruped_X(6:8, n), refQuadruped_Y(6:8, n), '.-', 'Color', [0 0.4470 0.7410], 'MarkerSize', 15, 'LineWidth', 1)% leg 2
    hold on
    plot(refQuadruped_X(9:11, n), refQuadruped_Y(9:11, n), '.-', 'Color', [0.8500 0.3250 0.0980], 'MarkerSize', 15, 'LineWidth', 2)% leg 3
    hold on
    plot(refQuadruped_X(12:14, n), refQuadruped_Y(12:14, n), '.-', 'Color', [0 0.4470 0.7410], 'MarkerSize', 15, 'LineWidth', 1)% leg 4
end
title('Ref Position')
axis equal
axisLimit_X = baseline_X;
axisLimit_Y = [-0.1; 1.1];
axis([axisLimit_X; axisLimit_Y]);

end

function [Quadruped_X, Quadruped_Y] = getQuadrupedConfiguration(plant, q)
x = q(1);
z = q(2);
xita_torso = q(3);
xita_thigh_1 = q(4);
xita_calf_1 = q(5);
xita_thigh_2 = q(6);
xita_calf_2 = q(7);
xita_thigh_3 = q(8);
xita_calf_3 = q(9);
xita_thigh_4 = q(10);
xita_calf_4 = q(11);

l_torso = plant.linkLength(1);
l_thigh = plant.linkLength(2);
l_calf = plant.linkLength(3);

% torso
torso_X = [x;...
    x + l_torso * sin(xita_torso)];
torso_Y = [z;...
    z - l_torso * cos(xita_torso)];

% leg 1
leg_1_X = [x; ...
    x + l_thigh * sin(xita_thigh_1);...
    x + l_thigh * sin(xita_thigh_1) + l_calf * sin(xita_calf_1)];
leg_1_Y = [z;...
    z - l_thigh * cos(xita_thigh_1);...
    z - l_thigh * cos(xita_thigh_1) - l_calf * cos(xita_calf_1)];

% leg 2
leg_2_X = [x;...
    x + l_thigh * sin(xita_thigh_2);...
    x + l_thigh * sin(xita_thigh_2) + l_calf * sin(xita_calf_2)];
leg_2_Y = [z;...
    z - l_thigh * cos(xita_thigh_2);...
    z - l_thigh * cos(xita_thigh_2) - l_calf * cos(xita_calf_2)];

% leg 3
leg_3_X = [x + l_torso * sin(xita_torso);...
    x + l_torso * sin(xita_torso) + l_thigh * sin(xita_thigh_3);...
    x + l_torso * sin(xita_torso) + l_thigh * sin(xita_thigh_3) + l_calf * sin(xita_calf_3)];
leg_3_Y = [z - l_torso * cos(xita_torso);...
    z - l_torso * cos(xita_torso) - l_thigh * cos(xita_thigh_3);...
    z - l_torso * cos(xita_torso) - l_thigh * cos(xita_thigh_3) - l_calf * cos(xita_calf_3)];

% leg 4
leg_4_X = [x + l_torso * sin(xita_torso);...
    x + l_torso * sin(xita_torso) + l_thigh * sin(xita_thigh_4);...
    x + l_torso * sin(xita_torso) + l_thigh * sin(xita_thigh_4) + l_calf * sin(xita_calf_4)];
leg_4_Y = [z - l_torso * cos(xita_torso);...
    z - l_torso * cos(xita_torso) - l_thigh * cos(xita_thigh_4);...
    z - l_torso * cos(xita_torso) - l_thigh * cos(xita_thigh_4) - l_calf * cos(xita_calf_4)];

Quadruped_X = [torso_X; leg_1_X; leg_2_X; leg_3_X; leg_4_X];
Quadruped_Y = [torso_Y; leg_1_Y; leg_2_Y; leg_3_Y; leg_4_Y];
end