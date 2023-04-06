function animateTrajectory(plant, simuTimeStep, InitState, tau, x, p)
%animateTrajectory
%   Detailed explanation goes here

DynVarLimit = plant.DynVarLimit;
Dim = plant.Dim;
% define time axis
nStages = size(x, 2);
timeAxis = 0 : simuTimeStep : nStages * simuTimeStep;

% get base line
baseline_X = [DynVarLimit.x_Min(1, 1); DynVarLimit.x_Max(1, 1)];  
baseline_Y = [0;0];

% get hopper position sequence based on given x sequence
Hopper_X = zeros(2, nStages + 1);
Hopper_Y = zeros(2, nStages + 1);

for n = 1 : nStages + 1
    if n == 1
        x_n = InitState;
    else
        x_n = x(:, n - 1);
    end
    [Hopper_X(:, n), Hopper_Y(:, n)] = getHopperConfiguration(x_n(1:4, :));
end

% compute value of friction and tangential velocity, and generate a trajectory-time sequence about contact force p for animation
K = zeros(Dim.p, nStages);
for n = 1 : nStages
    K(:, n) = plant.computeVI_Function(tau(:, n), x(:, n), p(:, n));
end

pN = [p(1, 1), p(1, :)];
pT = [tau(3, 1), tau(3, :)];
velT = [p(2, 1) - p(3, 1), p(2, :) - p(3, :)];

trajbody_X = cell(1, nStages + 1);
trajbody_Y = cell(1, nStages + 1);
trajfoot_X = cell(1, nStages + 1);
trajfoot_Y = cell(1, nStages + 1);

trajpN_X = cell(1, nStages + 1);
trajpN_Y = cell(1, nStages + 1);

trajpT_X = cell(1, nStages + 1);
trajpT_Y = cell(1, nStages + 1);
trajvelT_X = cell(1, nStages + 1);
trajvelT_Y = cell(1, nStages + 1);

for n = 1 : nStages + 1
    timeAxisSequence = [timeAxis(1:n), repmat(timeAxis(n), 1, nStages + 1 - n)];
    
    trajbody_X{1, n} = [Hopper_X(1, 1:n), repmat(Hopper_X(1, n), 1, nStages + 1 -n)];
    trajbody_Y{1, n} = [Hopper_Y(1, 1:n), repmat(Hopper_Y(1, n), 1, nStages + 1 -n)];
    
    trajfoot_X{1, n} = [Hopper_X(2, 1:n), repmat(Hopper_X(2, n), 1, nStages + 1 -n)];
    trajfoot_Y{1, n} = [Hopper_Y(2, 1:n), repmat(Hopper_Y(2, n), 1, nStages + 1 -n)];    
    
    trajpN_X{1, n} = timeAxisSequence;
    trajpN_Y{1, n} = [pN(1, 1:n), repmat(pN(1, n), 1, nStages + 1 -n)];
    
    trajpT_X{1, n} = timeAxisSequence;
    trajpT_Y{1, n} = [pT(1, 1:n), repmat(pT(1, n), 1, nStages + 1 - n)];
    
    trajvelT_X{1, n} = timeAxisSequence;
    trajvelT_Y{1, n} = [velT(1, 1:n), repmat(velT(1, n), 1, nStages + 1 - n)];    
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
% subplot(4,5,[1, 2, 3, 6, 7, 8, 11, 12, 13, 16, 17, 18])
subplot(1, 1, 1)
plot(baseline_X, baseline_Y, '.-k', 'MarkerSize', 1, 'LineWidth', 2);
hold on
% init, some middle and end state
MS_ID = [1; nStages + 1];
for i = 1 : size(MS_ID, 1)
    MSi = MS_ID(i);
    plot(Hopper_X(:, MSi), Hopper_Y(:, MSi), '.-', 'Color', [0 0.4470 0.7410], 'MarkerSize', 15, 'LineWidth', 1)
    hold on
end
% trajectory
Hopper = plot(Hopper_X(:, 1), Hopper_Y(:, 1), '.-', 'Color', [0 0.4470 0.7410], 'MarkerSize', 15, 'LineWidth', 1);
hold on
trajbody = plot(trajbody_X{1, 1}, trajbody_Y{1, 1}, '.--g', 'MarkerSize', 1, 'LineWidth', 1);
hold on
trajfoot = plot(trajfoot_X{1, 1}, trajfoot_Y{1, 1}, '.--r', 'MarkerSize', 1, 'LineWidth', 1);
hold on

% axis equal
axisLimit_X = baseline_X;
axisLimit_Y = [-0.1; 0.7];
% axisLimit_Y = [DynVarLimit.x_Min(2, 1) - DynVarLimit.x_Max(4, 1);...
%              DynVarLimit.x_Max(2, 1)]; 
axis([axisLimit_X; axisLimit_Y]);
xlabel('x_b [m]')
ylabel('y_b [m]')
timePrint = title(sprintf('Time: %0.2f sec', timeAxis(1)));

% subplot(4,5,[4, 5, 9, 10])
% trajpN = plot(trajpN_X{1, 1}, trajpN_Y{1, 1}, 'k', 'LineWidth', 1.2);
% axis([0; timeAxis(end); min(pN); max(pN)]);
% legend('pN', 'Location', 'northwest') 
% title('normal contact')

% subplot(4,5,[14, 15, 19, 20])
% trajpT = plot(trajpT_X{1, 1}, trajpT_Y{1, 1}, 'k', 'LineWidth', 1.2);
% hold on 
% trajvelT = plot(trajvelT_X{1, 1}, trajvelT_Y{1, 1}, 'b', 'LineWidth', 1.2);
% axis([0; timeAxis(end); min([pT,velT]); max([pT,velT])]);
% legend('pT', 'vel', 'Location', 'northwest') 
% xlabel('time(s)')
% title('friction contact')

%% animate trajectory
for n = 1 : nStages + 1
    % update XData and YData
    set(timePrint, 'String', sprintf('Time: %0.2f sec', timeAxis(n)));
    set(Hopper, 'XData', Hopper_X(:, n), 'YData', Hopper_Y(:, n));
    set(trajbody, 'XData', trajbody_X{1, n}, 'YData', trajbody_Y{1, n});
    set(trajfoot, 'XData', trajfoot_X{1, n}, 'YData', trajfoot_Y{1, n});
%     set(trajpN, 'XData', trajpN_X{1, n}, 'YData', trajpN_Y{1, n});    
%     set(trajpT, 'XData', trajpT_X{1, n}, 'YData', trajpT_Y{1, n})
%     set(trajvelT, 'XData', trajvelT_X{1, n}, 'YData', trajvelT_Y{1, n})
    
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
imwrite(mov, map, 'Hopper.gif', 'DelayTime', 0, 'LoopCount', inf)

end
%%
function [Hopper_X, Hopper_Y] = getHopperConfiguration(q)
Hopper_X = [q(1);...
    q(1) + q(4)*sin(q(3))];
Hopper_Y = [q(2);...
    q(2) - q(4)*cos(q(3))];
end
