function plotConfiguration(plant, InitState, RefState)
%plotConfiguration
%   Detailed explanation goes here

DynVarLimit = plant.DynVarLimit;
linkLength = plant.linkLength;
cartHeight = plant.cartHeight;
%
baseline_X = [DynVarLimit.x_Min(1, 1) - linkLength(1);...
              DynVarLimit.x_Max(1, 1) + linkLength(1)];
baseline_Y = [- 0.5 * cartHeight(1);...
              - 0.5 * cartHeight(1)];

[initCartPole_X, initCartPole_Y] = getCartPoleConfiguration(plant, InitState(1 : 2, 1));
[refCartPole_X, refCartPole_Y] = getCartPoleConfiguration(plant, RefState(1 : 2, 1));

figure(2000);
% plot InitPosition
subplot(2,1,1)
plot(baseline_X, baseline_Y, '.-k', 'MarkerSize', 1, 'LineWidth', 2);
hold on
plot(initCartPole_X(1 : 4, 1), initCartPole_Y(1 : 4, 1), '.-', 'Color', [0 0.4470 0.7410], 'MarkerSize', 2, 'LineWidth', 1);% cart
hold on
plot(initCartPole_X(5 : 6, 1), initCartPole_Y(5 : 6, 1), '.-', 'Color', [0.8500 0.3250 0.0980], 'MarkerSize', 2, 'LineWidth', 1);% pole
title('Init Position')
axis equal
axisLimit_X = baseline_X;
axisLimit_Y = [- 0.5 * cartHeight(1) - linkLength(1);...
               0.5 * cartHeight(1) + linkLength(1)]; 
axis([axisLimit_X; axisLimit_Y]);
% plot RefPosition
subplot(2,1,2)
plot(baseline_X, baseline_Y, '.-k', 'MarkerSize', 1, 'LineWidth', 2);
hold on
plot(refCartPole_X(1 : 4, 1), refCartPole_Y(1 : 4, 1), '.-', 'Color', [0 0.4470 0.7410],'MarkerSize', 2, 'LineWidth', 1);% cart
hold on
plot(refCartPole_X(5 : 6, 1), refCartPole_Y(5 : 6, 1), '.-', 'Color', [0.8500 0.3250 0.0980],'MarkerSize', 2, 'LineWidth', 1);% pole
title('Ref Position')
axis equal
axisLimit_X = baseline_X;
axisLimit_Y = [- 0.5 * cartHeight(1) - linkLength(1);...
               0.5 * cartHeight(1) + linkLength(1)]; 
axis([axisLimit_X; axisLimit_Y]);

end

%%
function [cartPole_X, cartPole_Y] = getCartPoleConfiguration(plant, q)
linkLength = plant.linkLength;
cartLength = plant.cartLength;
cartHeight = plant.cartHeight;
% cart coordinate
cart_X = [q(1) - 0.5 * cartLength(1);...
         q(1) - 0.5 * cartLength(1);...
         q(1) + 0.5 * cartLength(1);...
         q(1) + 0.5 * cartLength(1)];
cart_Y = [- 0.5 * cartHeight(1);...
          0.5 * cartHeight(1);...
          0.5 * cartHeight(1);...
         - 0.5 * cartHeight(1)];  
% pole coordinate
pole_X = [q(1);...
          q(1) + linkLength(1) * sin(q(2))];
pole_Y = [0;...
          - linkLength(1) * cos(q(2))];     
% cart-pole coordinate
cartPole_X = [cart_X;...
              pole_X];
cartPole_Y = [cart_Y;...
              pole_Y]; 

end
