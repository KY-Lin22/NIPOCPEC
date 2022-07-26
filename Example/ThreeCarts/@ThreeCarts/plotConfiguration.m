function plotConfiguration(plant, InitState, RefState)
%plotConfiguration
%   Detailed explanation goes here

DynVarLimit = plant.DynVarLimit;
cartLength = plant.cartLength;
%
baseline_X = [DynVarLimit.x_Min(1) - 0.5 * cartLength(1) - 1;...
              DynVarLimit.x_Max(3) + 0.5 * cartLength(3) + 1];
baseline_Y = [0; 0];

[initCart_X, initCart_Y] = getCartConfiguration(plant, InitState(1 : 3, 1));
[refCart_X, refCart_Y] = getCartConfiguration(plant, RefState(1 : 3, 1));

figure(2000);
subplot(2,1,1)
% plot InitPosition
plot(baseline_X, baseline_Y, '.-k', 'MarkerSize', 1, 'LineWidth', 2);
hold on
plot(initCart_X(1 : 4, 1), initCart_Y(1 : 4, 1), '.-r', 'MarkerSize', 2, 'LineWidth', 1);% cart 1
hold on
plot(initCart_X(5 : 8, 1), initCart_Y(5 : 8, 1), '.-g', 'MarkerSize', 2, 'LineWidth', 1);% cart 2
hold on
plot(initCart_X(9 : 12, 1), initCart_Y(9 : 12, 1), '.-b', 'MarkerSize', 2, 'LineWidth', 1);% cart 3
title('Init Position')
axis equal
axisLimit_X = baseline_X;
axisLimit_Y = [-1; 2 * plant.cartHeight(1)];
axis([axisLimit_X; axisLimit_Y]);
% plot RefPosition
subplot(2,1,2)
plot(baseline_X, baseline_Y, '.-k', 'MarkerSize', 1, 'LineWidth', 2);
hold on
plot(refCart_X(1 : 4, 1), refCart_Y(1 : 4, 1), '.-r','MarkerSize', 2, 'LineWidth', 1);% cart 1
hold on
plot(refCart_X(5 : 8, 1), refCart_Y(5 : 8, 1), '.-g',  'MarkerSize', 2, 'LineWidth', 1);% cart 2
hold on
plot(refCart_X(9 : 12, 1), refCart_Y(9 : 12, 1), '.-b', 'MarkerSize', 2, 'LineWidth', 1);% cart 3
title('Ref Position')
axis equal
axisLimit_X = baseline_X;
axisLimit_Y = [-1; 2 * plant.cartHeight(1)];
axis([axisLimit_X; axisLimit_Y]);

end

%%
function [cart_X, cart_Y] = getCartConfiguration(plant, q)
cartLength = plant.cartLength;
cartHeight = plant.cartHeight;
% cart 1 coordinate
cart_1_X = [q(1) - 0.5 * cartLength(1);...
            q(1) - 0.5 * cartLength(1);...
            q(1) + 0.5 * cartLength(1);...
            q(1) + 0.5 * cartLength(1)];
cart_1_Y = [0;...
            cartHeight(1);...
            cartHeight(1);...
            0];
% cart 2 coordinate
cart_2_X = [q(2) - 0.5 * cartLength(2);...
            q(2) - 0.5 * cartLength(2);...
            q(2) + 0.5 * cartLength(2);...
            q(2) + 0.5 * cartLength(2)];
cart_2_Y = [0;...
            cartHeight(2);...
            cartHeight(2);...
            0];    
% cart 3 coordinate        
cart_3_X = [q(3) - 0.5 * cartLength(3);...
            q(3) - 0.5 * cartLength(3);...
            q(3) + 0.5 * cartLength(3);...
            q(3) + 0.5 * cartLength(3)];
cart_3_Y = [0;...
            cartHeight(3);...
            cartHeight(3);...
            0];        
% cart coordinate
cart_X = [cart_1_X;...
          cart_2_X;...
          cart_3_X];
cart_Y = [cart_1_Y;...
          cart_2_Y;...
          cart_3_Y];      
end