function animateTraj_One_Oscillator(InitState, InitPhase, state, phase, nStages, timeStep)

%
state_extend = [InitState, state];
phase_extend = [InitPhase, phase];

timeAxis = 0 : timeStep : nStages * timeStep;

%
traj_X = cell(1, nStages + 1);
traj_Y = cell(1, nStages + 1);

timeAxisSequence = cell(1, nStages + 1);
trajp0_Y = cell(1, nStages + 1);
trajp1_Y = cell(1, nStages + 1);
trajphase_Y = cell(1, nStages + 1);
trajCosPhase_Y = cell(1, nStages + 1);

for n = 1 : nStages + 1
    timeAxisSequence{1, n} = [timeAxis(1:n), repmat(timeAxis(n), 1, nStages + 1 - n)];
    
    traj_X{1, n} = [state_extend(1, 1 : n), repmat(state_extend(1, n), 1, nStages + 1 - n)];
    traj_Y{1, n} = [state_extend(2, 1 : n), repmat(state_extend(2, n), 1, nStages + 1 - n)];
    
    trajp0_Y{1, n} = [state_extend(1, 1 : n), repmat(state_extend(1, n), 1, nStages + 1 - n)];
    trajp1_Y{1, n} = [state_extend(2, 1 : n), repmat(state_extend(2, n), 1, nStages + 1 - n)]; 
    trajphase_Y{1, n} = [phase_extend(1, 1 : n), repmat(phase_extend(1, n), 1, nStages + 1 - n)];
    trajCosPhase_Y{1, n} = [cos(phase_extend(1, 1 : n)), repmat(cos(phase_extend(1, n)), 1, nStages + 1 - n)];
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
subplot(4,6,[1,2,3,7,8,9,13,14,15,19,20,21])
trajstate = plot(state_extend(1, 1), state_extend(2, 1), '*', 'Color', [1 0 0], 'MarkerSize', 6, 'LineWidth', 2);
hold on
traj = plot(traj_X{1, 1}, traj_Y{1, 1}, '-', 'Color', [1 0 0], 'MarkerSize', 1, 'LineWidth', 1);
axisLimit_X = [min(state_extend(1, :)) - 0.5; max(state_extend(1, :)) + 0.5];
axisLimit_Y = [min(state_extend(2, :)) - 0.5; max(state_extend(2, :)) + 0.5];
axis equal
axis([axisLimit_X; axisLimit_Y]);
time = title(sprintf('time: %0.2f sec', timeAxis(1)));
xlabel('p0')
ylabel('p1')

subplot(4,6,[4,5,6])
trajp0 = plot(timeAxisSequence{1, 1}, trajp0_Y{1, 1}, '.--', 'Color', [0 0.4470 0.7410], 'LineWidth', 1);
axis([0; timeAxis(end); axisLimit_Y]);
ylabel('p0')

subplot(4,6,[10,11,12])
trajp1 = plot(timeAxisSequence{1, 1}, trajp1_Y{1, 1}, '.--', 'Color', [0 0.4470 0.7410], 'LineWidth', 1);
axis([0; timeAxis(end); axisLimit_Y]);
ylabel('p1')

subplot(4,6, [16,17,18])
trajphase = plot(timeAxisSequence{1, 1}, trajphase_Y{1, 1}, '.--', 'Color', [0 0.4470 0.7410], 'LineWidth', 1);
axis([0; timeAxis(end); [min(phase_extend(1, :)) - 0.5; max(phase_extend(1, :)) + 0.5]]);
ylabel('phase')

subplot(4,6, [22,23,24])
trajCosPhase = plot(timeAxisSequence{1, 1}, trajCosPhase_Y{1, 1}, '.--', 'Color', [0 0.4470 0.7410], 'LineWidth', 1);
axis([0; timeAxis(end); [-1.2; 1.2]]);
ylabel('cosPhase')
xlabel('time(s)')

for n = 1 : nStages + 1
    % update XData and YData   
    set(trajstate, 'XData', state_extend(1, n), 'YData', state_extend(2, n));
    set(traj, 'XData', traj_X{1, n}, 'YData', traj_Y{1, n});
    set(trajp0, 'XData', timeAxisSequence{1, n}, 'YData', trajp0_Y{1, n});
    set(trajp1, 'XData', timeAxisSequence{1, n}, 'YData', trajp1_Y{1, n});
    set(trajphase, 'XData', timeAxisSequence{1, n}, 'YData', trajphase_Y{1, n})
    set(trajCosPhase, 'XData', timeAxisSequence{1, n}, 'YData', trajCosPhase_Y{1, n})
    set(time, 'String', sprintf('time: %0.2f sec', timeAxis(n)));
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

