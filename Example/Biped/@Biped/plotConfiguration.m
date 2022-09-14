function plotConfiguration(plant, InitState, RefState)
%plotConfiguration
%   Detailed explanation goes here

DynVarLimit = plant.DynVarLimit;
%
baseline_X = [DynVarLimit.x_Min(1, 1) - 0.1; DynVarLimit.x_Max(1, 1) + 0.1];  
baseline_Y = [0;0];

[initBiped_X, initBiped_Y] = getBipedConfiguration(plant, InitState(1:7, 1));

refStageNum = size(RefState, 2);
refBiped_X = zeros(8, refStageNum);
refBiped_Y = zeros(8, refStageNum);
for n = 1 : refStageNum
    [refBiped_X(:, n), refBiped_Y(:, n)] = getBipedConfiguration(plant, RefState(1:7, n));
end

figure(2000)
% plot init position
subplot(2, 1, 1)
plot(baseline_X, baseline_Y, '-k', 'MarkerSize', 1, 'LineWidth', 2);
hold on
plot(initBiped_X(1:2, 1), initBiped_Y(1:2, 1), '.-k', 'MarkerSize', 15, 'LineWidth', 1)% torso
hold on
plot(initBiped_X(3:5, 1), initBiped_Y(3:5, 1), '.-', 'Color', [0.8500 0.3250 0.0980], 'MarkerSize', 15, 'LineWidth', 1)% leg 1
hold on
plot(initBiped_X(6:8, 1), initBiped_Y(6:8, 1), '.-', 'Color', [0 0.4470 0.7410], 'MarkerSize', 15, 'LineWidth', 1)% leg 2
title('Init Position')

axisLimit_X = baseline_X;
axisLimit_Y = [-0.1; 1];
axis([axisLimit_X; axisLimit_Y]);
% plot RefPosition
subplot(2,1,2)
plot(baseline_X, baseline_Y, '-k', 'MarkerSize', 1, 'LineWidth', 2);
for n = 1 : refStageNum
    hold on
    plot(refBiped_X(1:2, n), refBiped_Y(1:2, n), '.-k', 'MarkerSize', 15, 'LineWidth', 1)% torso
    hold on
    plot(refBiped_X(3:5, n), refBiped_Y(3:5, n), '.-', 'Color', [0.8500 0.3250 0.0980], 'MarkerSize', 15, 'LineWidth', 1)% leg 1
    hold on
    plot(refBiped_X(6:8, n), refBiped_Y(6:8, n), '.-', 'Color', [0 0.4470 0.7410], 'MarkerSize', 15, 'LineWidth', 1)% leg 2
end

title('Ref Position')
axisLimit_X = baseline_X;
axisLimit_Y = [-0.1; 1.2];
axis([axisLimit_X; axisLimit_Y]);

end

%%
function [Biped_X, Biped_Y] = getBipedConfiguration(plant, q)
x = q(1);
z = q(2);
xita_torso = q(3);
xita_thigh_1 = q(4);
xita_calf_1 = q(5);
xita_thigh_2 = q(6);
xita_calf_2 = q(7);

l_torso = plant.linkLength(1);
l_thigh = plant.linkLength(2);
l_calf = plant.linkLength(3);
% torso
torso_X = [x;...
    x - l_torso * sin(xita_torso)];
torso_Y = [z; ...
    z + l_torso * cos(xita_torso)];
% leg_1
leg_1_X = [x;...
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

Biped_X = [torso_X;...
    leg_1_X;...
    leg_2_X];
Biped_Y = [torso_Y;...
    leg_1_Y;...
    leg_2_Y];
end