function animateTrajectory_Acrobot(OCPEC, NLP, z_Opt)
%UNTITLED22 Summary of this function goes here
%   Detailed explanation goes here

nStages = OCPEC.nStages;
timeStep = OCPEC.timeStep;

Z_Opt = reshape(z_Opt, NLP.Dim.z_Node(end), OCPEC.nStages);

X_Opt = Z_Opt(1 : NLP.Dim.z_Node(1), :);
U_Opt = Z_Opt(NLP.Dim.z_Node(1) + 1 : NLP.Dim.z_Node(2), :);
LAMBDA_Opt = Z_Opt(NLP.Dim.z_Node(2) + 1 : NLP.Dim.z_Node(3), :);

F_FuncObj_map = OCPEC.FuncObj.F.map(nStages);
F_Opt = full(F_FuncObj_map(X_Opt, U_Opt, LAMBDA_Opt));

timeAxis = 0 : timeStep : nStages * timeStep;

% link length
L1 = 1;
L2 = 1;

% get link configuration (q1, q2) sequence based on given x_Opt sequence and generate a trajectory-time sequence for animation
link_position = [OCPEC.x0(1 : 2), X_Opt(1 : 2, :)];
Acrobot_X = zeros(8, nStages + 1);
Acrobot_Y = zeros(8, nStages + 1);
for n = 1 : nStages + 1
    [Acrobot_X(:, n), Acrobot_Y(:, n)] = getAcRobotConfiguration(link_position(:, n));
end

%%
% get figure size
figure(100)
% pos = get(gcf, 'Position');
% width = pos(3);
% height = pos(4);
% % define movie record
% mov = zeros(height*3/2, width*3/2, 1, nStages + 1, 'uint8');

% pre allocate
Acrobot = plot(Acrobot_X(1 : 4, 1), Acrobot_Y(1 : 4, 1), '.-k', 'MarkerSize', 10, 'LineWidth', 1.2);
hold on
maxLimit = plot(Acrobot_X(5 : 6, 1), Acrobot_Y(5 : 6, 1), '--b', 'MarkerSize', 10, 'LineWidth', 0.5);
hold on  
minLimit = plot(Acrobot_X(7 : 8, 1), Acrobot_Y(7 : 8, 1), '--b', 'MarkerSize', 10, 'LineWidth', 0.5);
hold on
axis equal
axisLimit_X = [-1 * (L1 + L2);...
               1 * (L1 + L2)];
axisLimit_Y = [-1 * (L1 + L2);...
               1 * (L1 + L2)];
axis([axisLimit_X; axisLimit_Y]);
timePrint = title(sprintf('Time: %0.2f sec', timeAxis(1)));

%% animate trajectory
animation = VideoWriter('Acrobot.mp4', 'MPEG-4');
open(animation);
for n = 1 : nStages + 1
    % update XData and YData
    set(Acrobot, 'XData', Acrobot_X(1 : 4, n), 'YData', Acrobot_Y(1 : 4, n));
    set(maxLimit, 'XData', Acrobot_X(5 : 6, n), 'YData', Acrobot_Y(5 : 6, n));
    set(minLimit, 'XData', Acrobot_X(7 : 8, n), 'YData', Acrobot_Y(7 : 8, n));
    set(timePrint, 'String', sprintf('Time: %0.2f sec', timeAxis(n)));
    % set(trajp1, 'XData', trajp1_X{1, n}, 'YData', trajp1_Y{1, n});
    % set(trajp2, 'XData', trajp2_X{1, n}, 'YData', trajp2_Y{1, n});
    % get frame as an image
    frame = getframe(gcf);
    writeVideo(animation, frame);
end


end



%%
function [AcRobot_X, AcRobot_Y] = getAcRobotConfiguration(q)
L1 = 1;
L2 = 1;
q2_min = -1/2*pi; % joint 2 max and min angle
q2_max = 1/2*pi;
% link 1 coordinate
link_1_X = [0;...
            L1 * sin(q(1))];
link_1_Y = [0;...
            -L1 * cos(q(1))];
% link 2 coordinate
link_2_X = [L1 * sin(q(1));...
            L1 * sin(q(1)) + L2 * sin(q(1) + q(2))];
link_2_Y = [-L1 * cos(q(1));...
            -L1 * cos(q(1)) - L2 * cos(q(1) + q(2))];
%        
maxLimit_X = [L1 * sin(q(1));...
              L1 * sin(q(1)) + 0.5 * L2 * sin(q(1) + q2_max)];   
maxLimit_Y = [-L1 * cos(q(1));...
              -L1 * cos(q(1)) - 0.5 * L2 * cos(q(1) + q2_max)];          
          
minLimit_X = [L1 * sin(q(1));...
              L1 * sin(q(1)) + 0.5 * L2 * sin(q(1) + q2_min)];   
minLimit_Y = [-L1 * cos(q(1));...
              -L1 * cos(q(1)) - 0.5 * L2 * cos(q(1) + q2_min)]; 
% link coordinate
AcRobot_X = [link_1_X; link_2_X; maxLimit_X; minLimit_X];
AcRobot_Y = [link_1_Y; link_2_Y; maxLimit_Y; minLimit_Y];
end
