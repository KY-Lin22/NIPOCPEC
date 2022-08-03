function animateTrajectory(plant, simuTimeStep, InitState, tau, x, p)
%animateTrajectory
%   Detailed explanation goes here

DynVarLimit = plant.DynVarLimit;
linkLength = plant.linkLength;
cartHeight = plant.cartHeight;
Dim = plant.Dim;
% define time axis
nStages = size(x, 2);
timeAxis = 0 : simuTimeStep : nStages * simuTimeStep;

% get baseline coordinate
baseline_X = [DynVarLimit.x_Min(1, 1) - linkLength(1);...
              DynVarLimit.x_Max(1, 1) + linkLength(1)];
baseline_Y = [- 0.5 * cartHeight(1);...
              - 0.5 * cartHeight(1)];
% get cartPole configuration sequence based on given x sequence
cartPole_X = zeros(6, nStages + 1);
cartPole_Y = zeros(6, nStages + 1);
for n = 1 : nStages + 1
    if n == 1
        x_n = InitState;
    else
        x_n = x(:, n - 1);
    end
    [cartPole_X(:, n), cartPole_Y(:, n)] = getCartPoleConfigurationC(plant, x_n(1:2, :));
end

% compute value of friction and tangential velocity, and generate a trajectory-time sequence for animation
K = zeros(Dim.p, nStages);
for n = 1 : nStages
    K(:, n) = plant.computeVI_Function(tau(:, n), x(:, n), p(:, n));
end
p = [p(:, 1), p];
K = [K(:, 1), K];
tau = [tau(:, 1), tau];% extend
trajPoleEnd_X = cell(1, nStages + 1);
trajPoleEnd_Y = cell(1, nStages + 1);
trajtau_X = cell(1, nStages + 1);
trajtau_Y = cell(1, nStages + 1);
trajp_X = cell(1, nStages + 1);
trajp_Y = cell(1, nStages + 1);
trajK_X = cell(1, nStages + 1);
trajK_Y = cell(1, nStages + 1);
for n = 1 : nStages + 1
    trajPoleEnd_X{1, n} = [cartPole_X(6, 1:n), repmat(cartPole_X(6, n), 1, nStages + 1 - n)];
    trajPoleEnd_Y{1, n} = [cartPole_Y(6, 1:n), repmat(cartPole_Y(6, n), 1, nStages + 1 - n)];
    
    trajtau_X{1, n} = [timeAxis(1:n), repmat(timeAxis(n), 1, nStages + 1 - n)];
    trajtau_Y{1, n} = [tau(1, 1:n), repmat(tau(1, n), 1, nStages + 1 - n)];
    trajp_X{1, n} = [timeAxis(1:n), repmat(timeAxis(n), 1, nStages + 1 - n)];
    trajp_Y{1, n} = [p(1, 1:n), repmat(p(1, n), 1, nStages + 1 - n)];
    trajK_X{1, n} = [timeAxis(1:n), repmat(timeAxis(n), 1, nStages + 1 - n)];
    trajK_Y{1, n} = [K(1, 1:n), repmat(K(1, n), 1, nStages + 1 - n)];
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
subplot(1,1,1)
% base
plot(baseline_X, baseline_Y, '.-k', 'MarkerSize', 1, 'LineWidth', 2);
hold on
% init state
plot(cartPole_X(1 : 4, 1), cartPole_Y(1 : 4, 1), '.-', 'Color', [0 0.4470 0.7410], 'MarkerSize', 2, 'LineWidth', 1);
hold on
plot(cartPole_X(5 : 6, 1), cartPole_Y(5 : 6, 1), '.-', 'Color', [0.8500 0.3250 0.0980], 'MarkerSize', 2, 'LineWidth', 1);
hold on
% some middle state
% MS1 = 24;
% plot(cartPole_X(1 : 4, MS1), cartPole_Y(1 : 4, MS1), '.-', 'Color', [0 0.4470 0.7410], 'MarkerSize', 2, 'LineWidth', 1);
% hold on
% plot(cartPole_X(5 : 6, MS1), cartPole_Y(5 : 6, MS1), '.-', 'Color', [0.8500 0.3250 0.0980], 'MarkerSize', 2, 'LineWidth', 1);
% hold on
% %
% MS2 = 38;
% plot(cartPole_X(1 : 4, MS2), cartPole_Y(1 : 4, MS2), '.-', 'Color', [0 0.4470 0.7410], 'MarkerSize', 2, 'LineWidth', 1);
% hold on
% plot(cartPole_X(5 : 6, MS2), cartPole_Y(5 : 6, MS2), '.-', 'Color', [0.8500 0.3250 0.0980], 'MarkerSize', 2, 'LineWidth', 1);
% hold on
% %
% MS3 = 50;
% plot(cartPole_X(1 : 4, MS3), cartPole_Y(1 : 4, MS3), '.-', 'Color', [0 0.4470 0.7410], 'MarkerSize', 2, 'LineWidth', 1);
% hold on
% plot(cartPole_X(5 : 6, MS3), cartPole_Y(5 : 6, MS3), '.-', 'Color', [0.8500 0.3250 0.0980], 'MarkerSize', 2, 'LineWidth', 1);
% hold on

% end state
plot(cartPole_X(1 : 4, end), cartPole_Y(1 : 4, end), '.-', 'Color', [0 0.4470 0.7410], 'MarkerSize', 2, 'LineWidth', 1);
hold on
plot(cartPole_X(5 : 6, end), cartPole_Y(5 : 6, end), '.-', 'Color', [0.8500 0.3250 0.0980], 'MarkerSize', 2, 'LineWidth', 1);
hold on
% trajectory
cart = plot(cartPole_X(1 : 4, 1), cartPole_Y(1 : 4, 1), '.-', 'Color', [0 0.4470 0.7410], 'MarkerSize', 2, 'LineWidth', 1);
hold on
pole = plot(cartPole_X(5 : 6, 1), cartPole_Y(5 : 6, 1), '.-', 'Color', [0.8500 0.3250 0.0980], 'MarkerSize', 2, 'LineWidth', 1);
hold on
PoleEnd = plot(trajPoleEnd_X{1, 1}, trajPoleEnd_Y{1, 1}, '.--r', 'MarkerSize', 1, 'LineWidth', 1);
xlabel('x_c [m]')
ylabel('y_c [m]')
axis equal
axisLimit_X = [-1; 7];
axisLimit_Y = [-1.2; 1.2];
% axisLimit_Y = [- 0.5 * cartHeight(1) - linkLength(1);...
%                0.5 * cartHeight(1) + linkLength(1)]; 
axis([axisLimit_X; axisLimit_Y]);

timePrint = title(sprintf('Time: %0.2f sec', timeAxis(1)));
% subplot(3,1,2)
% trajtau = plot(trajtau_X{1, 1}, trajtau_Y{1, 1}, 'k', 'LineWidth', 1.2);
% axis([0; timeAxis(end); min(tau); max(tau)]);
% xlabel('time [s]')
% title('control')
% 
% subplot(3,1,3)
% trajp = plot(trajp_X{1, 1}, trajp_Y{1, 1}, 'k', 'LineWidth', 1.2);
% hold on
% trajK = plot(trajK_X{1, 1}, trajK_Y{1, 1}, 'b', 'LineWidth', 1.2);
% axis([0; timeAxis(end); min([p,K]); max([p,K])]);
% legend('friction', 'cart vel') 
% xlabel('time [s]')
% title('equilibrium dynamics')
%% animate trajectory
for n = 1 : nStages + 1
    % update XData and YData   
    set(cart, 'XData', cartPole_X(1 : 4, n), 'YData', cartPole_Y(1 : 4, n));
    set(pole, 'XData', cartPole_X(5 : 6, n), 'YData', cartPole_Y(5 : 6, n));
    set(PoleEnd, 'XData', trajPoleEnd_X{1, n}, 'YData', trajPoleEnd_Y{1, n});
    set(timePrint, 'String', sprintf('Time: %0.2f sec', timeAxis(n)));
%     set(trajtau, 'XData', trajtau_X{1, n}, 'YData', trajtau_Y{1, n});
%     set(trajp, 'XData', trajp_X{1, n}, 'YData', trajp_Y{1, n});
%     set(trajK, 'XData', trajK_X{1, n}, 'YData', trajK_Y{1, n}); 
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
imwrite(mov, map, 'CartPoleWithFriction.gif', 'DelayTime', 0, 'LoopCount', inf)
end

%%
function [cartPole_X, cartPole_Y] = getCartPoleConfigurationC(plant, q)
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