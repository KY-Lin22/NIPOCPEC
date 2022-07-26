function animateTrajectory(plant, simuTimeStep, InitState, tau, x, p)
%animateTrajectory
%   Detailed explanation goes here

DynVarLimit = plant.DynVarLimit;
cartLength = plant.cartLength;

% define time axis
nStages = size(x, 2);
timeAxis = 0 : simuTimeStep : nStages * simuTimeStep;

% get baseline coordinate
baseline_X = [DynVarLimit.x_Min(1) - 0.5 * cartLength(1) - 2;...
              DynVarLimit.x_Max(3) + 0.5 * cartLength(3) + 2];
baseline_Y = [0; 0];

% get cart position sequence based on given x sequence
cart_X = zeros(12, nStages + 1);
cart_Y = zeros(12, nStages + 1);
for n = 1 : nStages + 1
    if n == 1
        x_n = InitState;
    else
        x_n = x(:, n - 1);
    end
    [cart_X(:, n), cart_Y(:, n)] = getCartConfiguration(plant, x_n(1 : 3, :));
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
subplot(2,1,1)
plot(baseline_X, baseline_Y, '.-k', 'MarkerSize', 1, 'LineWidth', 2);
hold on

cart_1 = plot(cart_X(1 : 4, 1), cart_Y(1 : 4, 1), '.-r', 'MarkerSize', 2, 'LineWidth', 1);
hold on
cart_2 = plot(cart_X(5 : 8, 1), cart_Y(5 : 8, 1), '.-g', 'MarkerSize', 2, 'LineWidth', 1);
hold on
cart_3 = plot(cart_X(9 : 12, 1), cart_Y(9 : 12, 1), '.-b', 'MarkerSize', 2, 'LineWidth', 1);
hold on

axisLimit_X = baseline_X;
axisLimit_Y = [-1; 3 * plant.cartHeight(1)];
axis([axisLimit_X; axisLimit_Y]);
timePrint = title(sprintf('Time: %0.2f sec', timeAxis(1)));

subplot(2,1,2)
trajp1 = plot(trajp1_X{1, 1}, trajp1_Y{1, 1},'r', 'LineWidth', 1.2);
hold on
trajp2 = plot(trajp2_X{1, 1}, trajp2_Y{1, 1}, 'b', 'LineWidth', 1.2);
axis([0; timeAxis(end); min([p(1, :), p(2, :)]); max([p(1, :), p(2, :)])]);
legend('cart 2 -> 1', 'cart 3 -> 2') 
xlabel('time(s)')
ylabel('force(N)')
title('normal contact force')

%% animate trajectory
for n = 1 : nStages + 1
    % update XData and YData   
    set(cart_1, 'XData', cart_X(1 : 4, n), 'YData', cart_Y(1 : 4, n));
    set(cart_2, 'XData', cart_X(5 : 8, n), 'YData', cart_Y(5 : 8, n));
    set(cart_3, 'XData', cart_X(9 : 12, n), 'YData', cart_Y(9 : 12, n));
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
imwrite(mov, map, 'ThreeCarts.gif', 'DelayTime', 0, 'LoopCount', inf)

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