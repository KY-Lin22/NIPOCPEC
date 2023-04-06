function animateTrajectory(plant, simuTimeStep, InitState, tau, x, p)
%animateTrajectory
%   Detailed explanation goes here

DynVarLimit = plant.DynVarLimit;

% define time axis
nStages = size(x, 2);
timeAxis = 0 : simuTimeStep : nStages * simuTimeStep;

% get base line
baseline_X = [DynVarLimit.x_Min(1, 1) - 0.2; DynVarLimit.x_Max(1, 1) + 0.1];  
baseline_Y = [0;0];

% get Quadruped position sequence based on given x sequence
Quadruped_X = zeros(14, nStages + 1);
Quadruped_Y = zeros(14, nStages + 1);
for n = 1 : nStages + 1
    if n == 1
        x_n = InitState;
    else
        x_n = x(:, n-1);
    end
    [Quadruped_X(:, n), Quadruped_Y(:, n)] = getQuadrupedConfiguration(plant, x_n(1:11, 1));
end
% compute value of friction and tangential velocity, and generate a trajectory-time sequence about contact force p for animation
pN_1 = [p(1, 1), p(1, :)];
pT_1 = [tau(9, 1), tau(9, :)];

pN_2 = [p(4, 1), p(4, :)];
pT_2 = [tau(10, 1), tau(10, :)];

pN_3 = [p(7, 1), p(7, :)];
pT_3 = [tau(11, 1), tau(11, :)];

pN_4 = [p(10, 1), p(10, :)];
pT_4 = [tau(12, 1), tau(12, :)];

trajbody_X = cell(1, nStages + 1);
trajbody_Y = cell(1, nStages + 1);

trajfoot_1_X = cell(1, nStages + 1);
trajfoot_1_Y = cell(1, nStages + 1);
trajfoot_2_X = cell(1, nStages + 1);
trajfoot_2_Y = cell(1, nStages + 1);

trajfoot_3_X = cell(1, nStages + 1);
trajfoot_3_Y = cell(1, nStages + 1);
trajfoot_4_X = cell(1, nStages + 1);
trajfoot_4_Y = cell(1, nStages + 1);

trajpN_1_X = cell(1, nStages + 1);
trajpN_1_Y = cell(1, nStages + 1);
trajpT_1_X = cell(1, nStages + 1);
trajpT_1_Y = cell(1, nStages + 1);

trajpN_2_X = cell(1, nStages + 1);
trajpN_2_Y = cell(1, nStages + 1);
trajpT_2_X = cell(1, nStages + 1);
trajpT_2_Y = cell(1, nStages + 1);

trajpN_3_X = cell(1, nStages + 1);
trajpN_3_Y = cell(1, nStages + 1);
trajpT_3_X = cell(1, nStages + 1);
trajpT_3_Y = cell(1, nStages + 1);

trajpN_4_X = cell(1, nStages + 1);
trajpN_4_Y = cell(1, nStages + 1);
trajpT_4_X = cell(1, nStages + 1);
trajpT_4_Y = cell(1, nStages + 1);

for n = 1 : nStages + 1
    timeAxisSequence = [timeAxis(1:n), repmat(timeAxis(n), 1, nStages + 1 - n)];
    trajbody_X{1, n} = [Quadruped_X(1, 1:n), repmat(Quadruped_X(1, n), 1, nStages + 1 -n)];
    trajbody_Y{1, n} = [Quadruped_Y(1, 1:n), repmat(Quadruped_Y(1, n), 1, nStages + 1 -n)];   
    
    trajfoot_1_X{1, n} = [Quadruped_X(5, 1:n), repmat(Quadruped_X(5, n), 1, nStages + 1 -n)];
    trajfoot_1_Y{1, n} = [Quadruped_Y(5, 1:n), repmat(Quadruped_Y(5, n), 1, nStages + 1 -n)];    
    trajfoot_2_X{1, n} = [Quadruped_X(8, 1:n), repmat(Quadruped_X(8, n), 1, nStages + 1 -n)];
    trajfoot_2_Y{1, n} = [Quadruped_Y(8, 1:n), repmat(Quadruped_Y(8, n), 1, nStages + 1 -n)]; 
    trajfoot_3_X{1, n} = [Quadruped_X(11, 1:n), repmat(Quadruped_X(11, n), 1, nStages + 1 -n)];
    trajfoot_3_Y{1, n} = [Quadruped_Y(11, 1:n), repmat(Quadruped_Y(11, n), 1, nStages + 1 -n)];    
    trajfoot_4_X{1, n} = [Quadruped_X(14, 1:n), repmat(Quadruped_X(14, n), 1, nStages + 1 -n)];
    trajfoot_4_Y{1, n} = [Quadruped_Y(14, 1:n), repmat(Quadruped_Y(14, n), 1, nStages + 1 -n)]; 
    
    trajpN_1_X{1, n} = timeAxisSequence;
    trajpN_1_Y{1, n} = [pN_1(1, 1:n), repmat(pN_1(1, n), 1, nStages + 1 -n)];
    trajpT_1_X{1, n} = timeAxisSequence;
    trajpT_1_Y{1, n} = [pT_1(1, 1:n), repmat(pT_1(1, n), 1, nStages + 1 - n)];    
    trajpN_2_X{1, n} = timeAxisSequence;
    trajpN_2_Y{1, n} = [pN_2(1, 1:n), repmat(pN_2(1, n), 1, nStages + 1 -n)];
    trajpT_2_X{1, n} = timeAxisSequence;
    trajpT_2_Y{1, n} = [pT_2(1, 1:n), repmat(pT_2(1, n), 1, nStages + 1 - n)];     

    trajpN_3_X{1, n} = timeAxisSequence;
    trajpN_3_Y{1, n} = [pN_3(1, 1:n), repmat(pN_3(1, n), 1, nStages + 1 -n)];
    trajpT_3_X{1, n} = timeAxisSequence;
    trajpT_3_Y{1, n} = [pT_3(1, 1:n), repmat(pT_3(1, n), 1, nStages + 1 - n)];    
    trajpN_4_X{1, n} = timeAxisSequence;
    trajpN_4_Y{1, n} = [pN_4(1, 1:n), repmat(pN_4(1, n), 1, nStages + 1 -n)];
    trajpT_4_X{1, n} = timeAxisSequence;
    trajpT_4_Y{1, n} = [pT_4(1, 1:n), repmat(pT_4(1, n), 1, nStages + 1 - n)]; 
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
subplot(4,8,[1, 2, 3, 4, 9, 10, 11, 12, 17, 18, 19, 20, 25, 26, 27, 28])
base = plot(baseline_X, baseline_Y, '.-k', 'MarkerSize', 1, 'LineWidth', 2);
hold on
Quadruped_torso = plot(Quadruped_X(1:2, 1), Quadruped_Y(1:2, 1), '.-k', 'MarkerSize', 15, 'LineWidth', 6);% torso
hold on
Quadruped_leg1 = plot(Quadruped_X(3:5, 1), Quadruped_Y(3:5, 1), '.-', 'Color', [0.8500 0.3250 0.0980], 'MarkerSize', 15, 'LineWidth', 1);% leg 1
hold on
Quadruped_leg2 = plot(Quadruped_X(6:8, 1), Quadruped_Y(6:8, 1), '.-', 'Color', [0 0.4470 0.7410], 'MarkerSize', 15, 'LineWidth', 1);% leg 2
hold on
Quadruped_leg3 = plot(Quadruped_X(9:11, 1), Quadruped_Y(9:11, 1), '.-', 'Color', [0.8500 0.3250 0.0980], 'MarkerSize', 15, 'LineWidth', 1);% leg 3
hold on
Quadruped_leg4 = plot(Quadruped_X(12:14, 1), Quadruped_Y(12:14, 1), '.-', 'Color', [0 0.4470 0.7410], 'MarkerSize', 15, 'LineWidth', 1);% leg 4
hold on
trajbody = plot(trajbody_X{1, 1}, trajbody_Y{1, 1}, '.--g', 'MarkerSize', 1, 'LineWidth', 1);
hold on
trajfoot_1 = plot(trajfoot_1_X{1, 1}, trajfoot_1_Y{1, 1}, '.--', 'Color', [0.8500 0.3250 0.0980], 'MarkerSize', 1, 'LineWidth', 1);
hold on
trajfoot_2 = plot(trajfoot_2_X{1, 1}, trajfoot_2_Y{1, 1}, '.--', 'Color', [0 0.4470 0.7410],  'MarkerSize', 1, 'LineWidth', 1);
hold on
trajfoot_3 = plot(trajfoot_3_X{1, 1}, trajfoot_3_Y{1, 1}, '.--', 'Color', [0.8500 0.3250 0.0980], 'MarkerSize', 1, 'LineWidth', 1);
hold on
trajfoot_4 = plot(trajfoot_4_X{1, 1}, trajfoot_4_Y{1, 1}, '.--', 'Color', [0 0.4470 0.7410],  'MarkerSize', 1, 'LineWidth', 1);
axis equal
axisLimit_X = baseline_X;
axisLimit_Y = [-0.1; 1.5];
axis([axisLimit_X; axisLimit_Y]);
timePrint = title(sprintf('Time: %0.2f sec', timeAxis(1)));

Y_Axis_min = 1.5 * min([pN_1, pT_1, pN_2, pT_2, pN_3, pT_3, pN_4, pT_4]);
Y_Axis_max = 1.5 * max([pN_1, pT_1, pN_2, pT_2, pN_3, pT_3, pN_4, pT_4]);

subplot(4,8,[5, 6, 13, 14])
trajpN_4 = plot(trajpN_4_X{1, 1}, trajpN_4_Y{1, 1}, '.-', 'Color', [0 0.4470 0.7410], 'LineWidth', 0.8);
hold on
trajpT_4 = plot(trajpT_4_X{1, 1}, trajpT_4_Y{1, 1}, '.--', 'Color', [0 0.4470 0.7410], 'LineWidth', 0.8);
axis([0; timeAxis(end); Y_Axis_min; Y_Axis_max]);
legend('pN(leg4)', 'pT(leg4)') 
set(gca,'YAxisLocation', 'right') 

subplot(4,8,[7, 8, 15, 16])
trajpN_3 = plot(trajpN_3_X{1, 1}, trajpN_3_Y{1, 1}, '.-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 0.8);
hold on
trajpT_3 = plot(trajpT_3_X{1, 1}, trajpT_3_Y{1, 1}, '.--', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 0.8);
axis([0; timeAxis(end); Y_Axis_min; Y_Axis_max]);
legend('pN(leg3)', 'pT(leg3)') 
set(gca,'Yticklabel',[],'YAxisLocation', 'right') 


subplot(4,8,[21, 22, 29, 30])
trajpN_2 = plot(trajpN_2_X{1, 1}, trajpN_2_Y{1, 1}, '.-', 'Color', [0 0.4470 0.7410], 'LineWidth', 0.8);
hold on
trajpT_2 = plot(trajpT_2_X{1, 1}, trajpT_2_Y{1, 1}, '.--', 'Color', [0 0.4470 0.7410], 'LineWidth', 0.8);
axis([0; timeAxis(end); Y_Axis_min; Y_Axis_max]);
legend('pN(leg2)', 'pT(leg2)') 
set(gca,'YAxisLocation', 'right') 
xlabel('time(s)')

subplot(4,8,[23, 24, 31, 32])
trajpN_1 = plot(trajpN_1_X{1, 1}, trajpN_1_Y{1, 1}, '.-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 0.8);
hold on
trajpT_1 = plot(trajpT_1_X{1, 1}, trajpT_1_Y{1, 1}, '.--', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 0.8);
axis([0; timeAxis(end); Y_Axis_min; Y_Axis_max]);
legend('pN(leg1)', 'pT(leg1)') 
set(gca,'Yticklabel',[], 'YAxisLocation', 'right') 
xlabel('time(s)')
%% animate trajectory
for n = 1 : nStages + 1
    % update XData and YData
    set(timePrint, 'String', sprintf('Time: %0.2f sec', timeAxis(n)));
    set(base, 'XData', baseline_X, 'YData', baseline_Y);    
    set(Quadruped_torso, 'XData', Quadruped_X(1:2, n), 'YData', Quadruped_Y(1:2, n));
    set(Quadruped_leg1, 'XData', Quadruped_X(3:5, n), 'YData', Quadruped_Y(3:5, n));
    set(Quadruped_leg2, 'XData', Quadruped_X(6:8, n), 'YData', Quadruped_Y(6:8, n));    
    set(Quadruped_leg3, 'XData', Quadruped_X(9:11, n), 'YData', Quadruped_Y(9:11, n));
    set(Quadruped_leg4, 'XData', Quadruped_X(12:14, n), 'YData', Quadruped_Y(12:14, n));     
    
    set(trajbody, 'XData', trajbody_X{1, n}, 'YData', trajbody_Y{1, n});
    set(trajfoot_1, 'XData', trajfoot_1_X{1, n}, 'YData', trajfoot_1_Y{1, n});
    set(trajfoot_2, 'XData', trajfoot_2_X{1, n}, 'YData', trajfoot_2_Y{1, n});
    set(trajfoot_3, 'XData', trajfoot_3_X{1, n}, 'YData', trajfoot_3_Y{1, n});
    set(trajfoot_4, 'XData', trajfoot_4_X{1, n}, 'YData', trajfoot_4_Y{1, n});
    
    set(trajpN_1, 'XData', trajpN_1_X{1, n}, 'YData', trajpN_1_Y{1, n});    
    set(trajpT_1, 'XData', trajpT_1_X{1, n}, 'YData', trajpT_1_Y{1, n});
    set(trajpN_2, 'XData', trajpN_2_X{1, n}, 'YData', trajpN_2_Y{1, n});    
    set(trajpT_2, 'XData', trajpT_2_X{1, n}, 'YData', trajpT_2_Y{1, n});    
    
    set(trajpN_3, 'XData', trajpN_3_X{1, n}, 'YData', trajpN_3_Y{1, n});    
    set(trajpT_3, 'XData', trajpT_3_X{1, n}, 'YData', trajpT_3_Y{1, n});
    set(trajpN_4, 'XData', trajpN_4_X{1, n}, 'YData', trajpN_4_Y{1, n});    
    set(trajpT_4, 'XData', trajpT_4_X{1, n}, 'YData', trajpT_4_Y{1, n});  
    
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
imwrite(mov, map, 'Quadruped.gif', 'DelayTime', 0, 'LoopCount', inf)

end

%%
function [Quadruped_X, Quadruped_Y] = getQuadrupedConfiguration(plant, q)
x = q(1);
z = q(2);
xita_torso = q(3);
xita_thigh_1 = q(4);
xita_calf_1 = q(5);
xita_thigh_2 = q(6);
xita_calf_2 = q(7);
xita_thigh_3 = q(8);
xita_calf_3 = q(9);
xita_thigh_4 = q(10);
xita_calf_4 = q(11);

l_torso = plant.linkLength(1);
l_thigh = plant.linkLength(2);
l_calf = plant.linkLength(3);

% torso
torso_X = [x;...
    x + l_torso * sin(xita_torso)];
torso_Y = [z;...
    z - l_torso * cos(xita_torso)];

% leg 1
leg_1_X = [x; ...
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

% leg 3
leg_3_X = [x + l_torso * sin(xita_torso);...
    x + l_torso * sin(xita_torso) + l_thigh * sin(xita_thigh_3);...
    x + l_torso * sin(xita_torso) + l_thigh * sin(xita_thigh_3) + l_calf * sin(xita_calf_3)];
leg_3_Y = [z - l_torso * cos(xita_torso);...
    z - l_torso * cos(xita_torso) - l_thigh * cos(xita_thigh_3);...
    z - l_torso * cos(xita_torso) - l_thigh * cos(xita_thigh_3) - l_calf * cos(xita_calf_3)];

% leg 4
leg_4_X = [x + l_torso * sin(xita_torso);...
    x + l_torso * sin(xita_torso) + l_thigh * sin(xita_thigh_4);...
    x + l_torso * sin(xita_torso) + l_thigh * sin(xita_thigh_4) + l_calf * sin(xita_calf_4)];
leg_4_Y = [z - l_torso * cos(xita_torso);...
    z - l_torso * cos(xita_torso) - l_thigh * cos(xita_thigh_4);...
    z - l_torso * cos(xita_torso) - l_thigh * cos(xita_thigh_4) - l_calf * cos(xita_calf_4)];

Quadruped_X = [torso_X; leg_1_X; leg_2_X; leg_3_X; leg_4_X];
Quadruped_Y = [torso_Y; leg_1_Y; leg_2_Y; leg_3_Y; leg_4_Y];
end

