function animateTrajectory(plant, simuTimeStep, InitState, tau, x, p)
%animateTrajectory
%   Detailed explanation goes here

%%
DynVarLimit = plant.DynVarLimit;
Dim = plant.Dim;

% define time axis
nStages = size(x, 2);
timeAxis = 0 : simuTimeStep : nStages * simuTimeStep;

% get baseline coordinate
baseLine_horizontal_X = [DynVarLimit.x_Min(1); DynVarLimit.x_Max(1)];
baseLine_horizontal_Y = [0; 0];  
baseLine_vertical_X = [0; 0];
baseLine_vertical_Y = [DynVarLimit.x_Min(2); DynVarLimit.x_Max(2)]; 

% get state sequence based on given x sequence
state_X = [InitState(1), x(1, :)];
state_Y = [InitState(2), x(2, :)];

% generate a trajectory-time sequence for animation
traj_X = cell(1, nStages + 1);
traj_Y = cell(1, nStages + 1);
for n = 1 : nStages + 1
    traj_X{1, n} = [state_X(1, 1 : n), repmat(state_X(1, n), 1, nStages + 1 - n)];
    traj_Y{1, n} = [state_Y(1, 1 : n), repmat(state_Y(1, n), 1, nStages + 1 - n)];
end

% compute value of p, K and generate a trajectory-time sequence for animation
K = zeros(Dim.p, nStages);
for n = 1 : nStages
    K(:, n) = plant.computeVI_Function(tau(:, n), x(:, n), p(:, n));
end

p = [p(:, 1), p];
K = [K(:, 1), K];% extend
trajp_X = cell(1, nStages + 1);
trajp_Y = cell(1, nStages + 1);
trajK_X = cell(1, nStages + 1);
trajK_Y = cell(1, nStages + 1);
for n = 1 : nStages + 1
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
plot(baseLine_horizontal_X, baseLine_horizontal_Y, '.-k', 'MarkerSize', 0.5, 'LineWidth', 0.5);
hold on
plot(baseLine_vertical_X, baseLine_vertical_Y, '.-k', 'MarkerSize', 0.5, 'LineWidth', 0.5);
hold on
state = plot(state_X(1, 1), state_Y(1, 1), '*', 'Color', [1 0 0], 'MarkerSize', 6, 'LineWidth', 2);
hold on
traj = plot(traj_X{1, 1}, traj_Y{1, 1}, '-', 'Color', [1 0 0], 'MarkerSize', 1, 'LineWidth', 1);
axisLimit_X = baseLine_horizontal_X;
axisLimit_Y = baseLine_vertical_Y; 
axis([axisLimit_X; axisLimit_Y]);
timePrint = title(sprintf('Time: %0.2f sec', timeAxis(1)));
xlabel('x1')
ylabel('x2')

% subplot(2,1,2)
% trajp = plot(trajp_X{1, 1}, trajp_Y{1, 1}, 'k', 'LineWidth', 1.2);
% hold on
% trajK = plot(trajK_X{1, 1}, trajK_Y{1, 1}, 'b', 'LineWidth', 1.2);
% axis([0; timeAxis(end); min([p,K]); max([p,K])]);
% legend('p', 'K')
% xlabel('time(s)')
% title('equilibrium dynamics');

%% animate trajectory
for n = 1 : nStages + 1
    % update XData and YData   
    set(state, 'XData', state_X(1, n), 'YData', state_Y(1, n));
    set(traj, 'XData', traj_X{1, n}, 'YData', traj_Y{1, n});
    set(timePrint, 'String', sprintf('Time: %0.2f sec', timeAxis(n)));
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
imwrite(mov, map, 'AffineDVI.gif', 'DelayTime', 0, 'LoopCount', inf)

end

