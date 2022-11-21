function plotConfiguration(plant, InitState, RefState)
%plotConfiguration
%   Detailed explanation goes here

DynVarLimit = plant.DynVarLimit;

% 
baseline_X = [DynVarLimit.x_Min(1, 1) - 1; DynVarLimit.x_Max(1, 1) + 1];  
baseline_Y = [0;0];

[initHopper_X, initHopper_Y] = getHopperConfiguration(InitState(1:4, 1));

refStageNum = size(RefState, 2);
refHopper_X = zeros(2, refStageNum);
refHopper_Y = zeros(2, refStageNum);
for n = 1 : refStageNum
    [refHopper_X(:, n), refHopper_Y(:, n)] = getHopperConfiguration(RefState(1:4, n));
end

figure(2000);
% plot InitPosition
subplot(2,1,1)
plot(baseline_X, baseline_Y, '.-k', 'MarkerSize', 1, 'LineWidth', 2);
hold on
plot(initHopper_X, initHopper_Y, '.-', 'Color', [0 0.4470 0.7410], 'MarkerSize', 15, 'LineWidth', 1);
title('Init Position')
axis equal
axisLimit_X = baseline_X;
axisLimit_Y = [DynVarLimit.x_Min(2, 1) - DynVarLimit.x_Max(4, 1);...
             DynVarLimit.x_Max(2, 1) + DynVarLimit.x_Max(4, 1)]; 
axis([axisLimit_X; axisLimit_Y]);
% plot RefPosition
subplot(2,1,2)
plot(baseline_X, baseline_Y, '.-k', 'MarkerSize', 1, 'LineWidth', 2);
for n = 1 : refStageNum
   hold on
   plot(refHopper_X(:, n), refHopper_Y(:, n), '.-', 'Color', [0 0.4470 0.7410], 'MarkerSize', 15, 'LineWidth', 1);
end

title('Ref Position')
axis equal
axisLimit_X = baseline_X;
axisLimit_Y = [DynVarLimit.x_Min(2, 1) - DynVarLimit.x_Max(4, 1);...
             DynVarLimit.x_Max(2, 1) + DynVarLimit.x_Max(4, 1)]; 
axis([axisLimit_X; axisLimit_Y]);

end

%%
function [Hopper_X, Hopper_Y] = getHopperConfiguration(q)
Hopper_X = [q(1);...
    q(1) + q(4)*sin(q(3))];
Hopper_Y = [q(2);...
    q(2) - q(4)*cos(q(3))];
end
