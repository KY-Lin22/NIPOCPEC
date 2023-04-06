function animateSystemSimulationResult(plant, InitState, state, nStages, timeStep)
%animateSystemSimulationResult
%   Detailed explanation goes here
num_Osci = length(plant.InitPhase);
timeAxis = 0 : timeStep : nStages * timeStep;

state_extend = [InitState, state];
%
traj_X = cell(num_Osci, nStages + 1);
traj_Y = cell(num_Osci, nStages + 1);
trajfootpos_X = cell(num_Osci, nStages + 1);
trajfootpos_Y = cell(num_Osci, nStages + 1);

for n = 1 : nStages + 1
    for i = 1 : num_Osci
        state_extend_i = state_extend(1 + (i - 1) * 4 : 4 + (i - 1) * 4, :);
        traj_X{i, n} = [state_extend_i(1, 1 : n), repmat(state_extend_i(1, n), 1, nStages + 1 - n)];
        traj_Y{i, n} = [state_extend_i(2, 1 : n), repmat(state_extend_i(2, n), 1, nStages + 1 - n)];
        
        trajfootpos_X{i, n} = [state_extend_i(4, 1 : n), repmat(state_extend_i(4, n), 1, nStages + 1 - n)];
        trajfootpos_Y{i, n} = [state_extend_i(3, 1 : n), repmat(state_extend_i(3, n), 1, nStages + 1 - n)];
    end
end

%%
% get figure size
figure(100)
pos = get(gcf, 'Position');
width = pos(3);
height = pos(4);
% define movie record
mov = zeros(height*3/2, width*3/2, 1, nStages + 1, 'uint8');

% axis limit
axisLimit_X_Osci = zeros(num_Osci, 2);
axisLimit_Y_Osci = zeros(num_Osci, 2);
axisLimit_X_FootTraj = zeros(num_Osci, 2);
axisLimit_Y_FootTraj = zeros(num_Osci, 2);
for i = 1 : num_Osci
    state_extend_i = state_extend(1 + (i - 1) * 4 : 4 + (i - 1) * 4, :);
    % Oscillator
    axisLimit_X_Osci(i, :) = [min(state_extend_i(1, :)) - 0.5, max(state_extend_i(1, :)) + 0.5];
    axisLimit_Y_Osci(i, :) = [min(state_extend_i(2, :)) - 0.5, max(state_extend_i(2, :)) + 0.5];    
    % Foot Trajectory
    axisLimit_X_FootTraj(i, :) = [min(state_extend_i(4, :)) - 0.1, max(state_extend_i(4, :)) + 0.1];
    axisLimit_Y_FootTraj(i, :) = [min(state_extend_i(3, :)) - 0.1, max(state_extend_i(3, :)) + 0.1];       
end
axisLimit_Osci = [min(axisLimit_X_Osci(:, 1)), max(axisLimit_X_Osci(:, 2)),...
    min(axisLimit_Y_Osci(:, 1)), max(axisLimit_Y_Osci(:, 2))];
axisLimit_FootTraj = [min(axisLimit_X_FootTraj(:, 1)), max(axisLimit_X_FootTraj(:, 2)),...
    min(axisLimit_Y_FootTraj(:, 1)), max(axisLimit_Y_FootTraj(:, 2))];

% pre allocate
trajstate = cell(num_Osci, 1);
traj = cell(num_Osci, 1);
trajfootpos = cell(num_Osci, 1);
for i = 1 : num_Osci
    state_extend_i = state_extend(1 + (i - 1) * 4 : 4 + (i - 1) * 4, :);
    subplot(num_Osci, 2, 1 + (i - 1) * 2)

    trajstate_i = plot(state_extend_i(1, 1), state_extend_i(2, 1), '*', 'Color', [1 0 0], 'MarkerSize', 6, 'LineWidth', 2);
    trajstate{i, 1} = trajstate_i;
    if i == 1
        time = title(sprintf('time: %0.2f sec', timeAxis(1)));
    end
    hold on
    traj_i = plot(traj_X{i, 1}, traj_Y{i, 1}, '-', 'Color', [1 0 0], 'MarkerSize', 1, 'LineWidth', 1);
    traj{i, 1} = traj_i;
%     axis equal
    axis(axisLimit_Osci);
    xlabel('\rho_1')
    ylabel('\rho_2')
    
    subplot(num_Osci, 2, 2 + (i - 1) * 2)
    trajfootpos_i = plot(trajfootpos_X{i, 1}, trajfootpos_Y{i, 1}, '.--', 'Color', [0 0.4470 0.7410], 'LineWidth', 1);
    trajfootpos{i, 1} = trajfootpos_i;
    hold on
    plot(axisLimit_FootTraj(1:2), [0; 0])
%     axis equal
    axis(axisLimit_FootTraj);
end

%%
for n = 1 : nStages + 1
    % update XData and YData   
    set(time, 'String', sprintf('time: %0.2f sec', timeAxis(n)));
    for i = 1 : num_Osci
        state_extend_i = state_extend(1 + (i - 1) * 4 : 4 + (i - 1) * 4, :);
        set(trajstate{i, 1}, 'XData', state_extend_i(1, n), 'YData', state_extend_i(2, n));
        set(traj{i, 1}, 'XData', traj_X{i, n}, 'YData', traj_Y{i, n});
        set(trajfootpos{i, 1}, 'XData', trajfootpos_X{i, n}, 'YData', trajfootpos_Y{i, n});
    end

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
imwrite(mov, map, 'holf_oscillator.gif', 'DelayTime', 0, 'LoopCount', inf)

end

