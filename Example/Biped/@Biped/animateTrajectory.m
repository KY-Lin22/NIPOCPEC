function animateTrajectory(plant, simuTimeStep, InitState, tau, x, p)
%animateTrajectory
%   Detailed explanation goes here

DynVarLimit = plant.DynVarLimit;

% define time axis
nStages = size(x, 2);
timeAxis = 0 : simuTimeStep : nStages * simuTimeStep;

% get base line
baseline_X = [DynVarLimit.x_Min(1, 1) - 0.1; DynVarLimit.x_Max(1, 1) + 0.1];  
baseline_Y = [0;0];

% get biped position sequence based on given x sequence
Biped_X = zeros(8, nStages + 1);
Biped_Y = zeros(8, nStages + 1);
for n = 1 : nStages + 1
    if n == 1
        x_n = InitState;
    else
        x_n = x(:, n - 1);
    end
    [Biped_X(:, n), Biped_Y(:, n)] = getBipedConfiguration(plant, x_n(1:7, 1));
end

% compute value of friction and tangential velocity, and generate a trajectory-time sequence about contact force p for animation
pN_1 = [p(1, 1), p(1, :)];
pT_1 = [tau(6, 1), tau(6, :)];
% velT_1 = [p(2, 1) - p(3, 1), p(2, :) - p(3, :)];

pN_2 = [p(4, 1), p(4, :)];
pT_2 = [tau(7, 1), tau(7, :)];
% velT_2 = [p(5, 1) - p(6, 1), p(5, :) - p(6, :)];

trajbody_X = cell(1, nStages + 1);
trajbody_Y = cell(1, nStages + 1);
trajfoot_1_X = cell(1, nStages + 1);
trajfoot_1_Y = cell(1, nStages + 1);
trajfoot_2_X = cell(1, nStages + 1);
trajfoot_2_Y = cell(1, nStages + 1);

trajpN_1_X = cell(1, nStages + 1);
trajpN_1_Y = cell(1, nStages + 1);
trajpT_1_X = cell(1, nStages + 1);
trajpT_1_Y = cell(1, nStages + 1);

trajpN_2_X = cell(1, nStages + 1);
trajpN_2_Y = cell(1, nStages + 1);
trajpT_2_X = cell(1, nStages + 1);
trajpT_2_Y = cell(1, nStages + 1);

for n = 1 : nStages + 1
    timeAxisSequence = [timeAxis(1:n), repmat(timeAxis(n), 1, nStages + 1 - n)];
    
    trajbody_X{1, n} = [Biped_X(1, 1:n), repmat(Biped_X(1, n), 1, nStages + 1 -n)];
    trajbody_Y{1, n} = [Biped_Y(1, 1:n), repmat(Biped_Y(1, n), 1, nStages + 1 -n)];      
    trajfoot_1_X{1, n} = [Biped_X(5, 1:n), repmat(Biped_X(5, n), 1, nStages + 1 -n)];
    trajfoot_1_Y{1, n} = [Biped_Y(5, 1:n), repmat(Biped_Y(5, n), 1, nStages + 1 -n)];    
    trajfoot_2_X{1, n} = [Biped_X(8, 1:n), repmat(Biped_X(8, n), 1, nStages + 1 -n)];
    trajfoot_2_Y{1, n} = [Biped_Y(8, 1:n), repmat(Biped_Y(8, n), 1, nStages + 1 -n)]; 
    
    trajpN_1_X{1, n} = timeAxisSequence;
    trajpN_1_Y{1, n} = [pN_1(1, 1:n), repmat(pN_1(1, n), 1, nStages + 1 -n)];
    trajpT_1_X{1, n} = timeAxisSequence;
    trajpT_1_Y{1, n} = [pT_1(1, 1:n), repmat(pT_1(1, n), 1, nStages + 1 - n)];    
    trajpN_2_X{1, n} = timeAxisSequence;
    trajpN_2_Y{1, n} = [pN_2(1, 1:n), repmat(pN_2(1, n), 1, nStages + 1 -n)];
    trajpT_2_X{1, n} = timeAxisSequence;
    trajpT_2_Y{1, n} = [pT_2(1, 1:n), repmat(pT_2(1, n), 1, nStages + 1 - n)];       
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
subplot(4,6,[1, 2, 3, 4, 7, 8, 9, 10, 13, 14, 15, 16, 19, 20, 21, 22])
plot(baseline_X, baseline_Y, '.-k', 'MarkerSize', 1, 'LineWidth', 2);
hold on
Biped_torso = plot(Biped_X(1:2, 1), Biped_Y(1:2, 1), '.-k', 'MarkerSize', 15, 'LineWidth', 1);% torso
hold on
Biped_leg1 = plot(Biped_X(3:5, 1), Biped_Y(3:5, 1), '.-', 'Color', [0.8500 0.3250 0.0980], 'MarkerSize', 15, 'LineWidth', 1);% leg 1
hold on
Biped_leg2 = plot(Biped_X(6:8, 1), Biped_Y(6:8, 1), '.-', 'Color', [0 0.4470 0.7410], 'MarkerSize', 15, 'LineWidth', 1);% leg 2
hold on
trajbody = plot(trajbody_X{1, 1}, trajbody_Y{1, 1}, '.--g', 'MarkerSize', 1, 'LineWidth', 1);
hold on
trajfoot_1 = plot(trajfoot_1_X{1, 1}, trajfoot_1_Y{1, 1}, '.--', 'Color', [0.8500 0.3250 0.0980], 'MarkerSize', 1, 'LineWidth', 1);
hold on
trajfoot_2 = plot(trajfoot_2_X{1, 1}, trajfoot_2_Y{1, 1}, '.--', 'Color', [0 0.4470 0.7410],  'MarkerSize', 1, 'LineWidth', 1);
axis equal
axisLimit_X = baseline_X;
axisLimit_Y = [-0.1; 1];
axis([axisLimit_X; axisLimit_Y]);
timePrint = title(sprintf('Time: %0.2f sec', timeAxis(1)));

Y_Axis_min = 1.5 * min([pN_1, pT_1, pN_2, pT_2]);
Y_Axis_max = 1.5 * max([pN_1, pT_1, pN_2, pT_2]);

subplot(4,6,[5, 6, 11, 12])
trajpN_1 = plot(trajpN_1_X{1, 1}, trajpN_1_Y{1, 1}, '.-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1);
hold on
trajpT_1 = plot(trajpT_1_X{1, 1}, trajpT_1_Y{1, 1}, '.--', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1);
axis([0; timeAxis(end); Y_Axis_min; Y_Axis_max]);
legend('pN(leg1)', 'pT(leg1)') 

subplot(4,6,[17, 18, 23, 24])
trajpN_2 = plot(trajpN_2_X{1, 1}, trajpN_2_Y{1, 1}, '.-', 'Color', [0 0.4470 0.7410], 'LineWidth', 1);
hold on
trajpT_2 = plot(trajpT_2_X{1, 1}, trajpT_2_Y{1, 1}, '.--', 'Color', [0 0.4470 0.7410], 'LineWidth', 1);
axis([0; timeAxis(end); Y_Axis_min; Y_Axis_max]);
legend('pN(leg2)', 'pT(leg2)') 
xlabel('time(s)')

%% animate trajectory
for n = 1 : nStages + 1
    % update XData and YData
    set(timePrint, 'String', sprintf('Time: %0.2f sec', timeAxis(n)));   
    set(Biped_torso, 'XData', Biped_X(1:2, n), 'YData', Biped_Y(1:2, n));
    set(Biped_leg1, 'XData', Biped_X(3:5, n), 'YData', Biped_Y(3:5, n));
    set(Biped_leg2, 'XData', Biped_X(6:8, n), 'YData', Biped_Y(6:8, n));
    
    set(trajbody, 'XData', trajbody_X{1, n}, 'YData', trajbody_Y{1, n});
    set(trajfoot_1, 'XData', trajfoot_1_X{1, n}, 'YData', trajfoot_1_Y{1, n});
    set(trajfoot_2, 'XData', trajfoot_2_X{1, n}, 'YData', trajfoot_2_Y{1, n});
    
    set(trajpN_1, 'XData', trajpN_1_X{1, n}, 'YData', trajpN_1_Y{1, n});    
    set(trajpT_1, 'XData', trajpT_1_X{1, n}, 'YData', trajpT_1_Y{1, n});
    set(trajpN_2, 'XData', trajpN_2_X{1, n}, 'YData', trajpN_2_Y{1, n});    
    set(trajpT_2, 'XData', trajpT_2_X{1, n}, 'YData', trajpT_2_Y{1, n});   
    
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
imwrite(mov, map, 'Biped.gif', 'DelayTime', 0, 'LoopCount', inf)

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