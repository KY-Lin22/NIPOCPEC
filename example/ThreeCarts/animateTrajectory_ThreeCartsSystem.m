function animateTrajectory_ThreeCartsSystem(OCPEC, NLP, z_Opt)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%%
nStages = OCPEC.nStages;
timeStep = OCPEC.timeStep;

Z_Opt = reshape(z_Opt, NLP.Dim.z_Node(end), OCPEC.nStages);

X_Opt = Z_Opt(1 : NLP.Dim.z_Node(1), :);
U_Opt = Z_Opt(NLP.Dim.z_Node(1) + 1 : NLP.Dim.z_Node(2), :);
LAMBDA_Opt = Z_Opt(NLP.Dim.z_Node(2) + 1 : NLP.Dim.z_Node(3), :);

F_FuncObj_map = OCPEC.FuncObj.F.map(nStages);
F_Opt = full(F_FuncObj_map(X_Opt, U_Opt, LAMBDA_Opt));

timeAxis = 0 : timeStep : nStages * timeStep;

% cart length
L1 = 2; 
L2 = 2;
L3 = 2;

% get baseline coordinate
baseline_X = [-8; 8];
baseline_Y = [0; 0];

% get reference coordinate
cart_ref = [-5; 0; 5];
cart_1_ref = nsidedpoly(4, 'Center', [cart_ref(1), L1/2], 'SideLength', L1);
cart_2_ref = nsidedpoly(4, 'Center', [cart_ref(2), L2/2], 'SideLength', L2);
cart_3_ref = nsidedpoly(4, 'Center', [cart_ref(3), L3/2], 'SideLength', L3);

% get cart configuration (q1, q2, q3) sequence based on given x_Opt sequence and generate a trajectory-time sequence for animation
cart_position = [OCPEC.x0(1 : 3), X_Opt(1 : 3, :)];

cart_1_position = cell(1, nStages + 1);
cart_2_position = cell(1, nStages + 1);
cart_3_position = cell(1, nStages + 1);
for n = 1 : nStages + 1
    cart_1_position{1, n} = nsidedpoly(4, 'Center', [cart_position(1, n), L1/2], 'SideLength', L1);
    cart_2_position{1, n} = nsidedpoly(4, 'Center', [cart_position(2, n), L2/2], 'SideLength', L2);
    cart_3_position{1, n} = nsidedpoly(4, 'Center', [cart_position(3, n), L3/2], 'SideLength', L3);
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
plot(baseline_X, baseline_Y, '.-k', 'MarkerSize', 1, 'LineWidth', 2);
hold on
plot(cart_1_ref,  'FaceAlpha', 0.01, 'LineStyle', '--')
hold on
plot(cart_2_ref,  'FaceAlpha', 0.01, 'LineStyle', '--')
hold on
plot(cart_3_ref,  'FaceAlpha', 0.01, 'LineStyle', '--')
hold on
cart_1 = plot(cart_1_position{1, 1});
hold on
cart_2 = plot(cart_2_position{1, 1});
hold on
cart_3 = plot(cart_3_position{1, 1});
hold on
% axisLimit_X = baseline_X;
% axisLimit_Y = [-1; inf];
% axis([axisLimit_X; axisLimit_Y]);
axis equal
axisLimit_X = baseline_X;
axisLimit_Y = [-1; 3];
axis([axisLimit_X; axisLimit_Y]);
timePrint = title(sprintf('Time: %0.2f sec', timeAxis(1)));

%% animate trajectory
animation = VideoWriter('ThreeCarts.mp4', 'MPEG-4');
open(animation);
for n = 1 : nStages + 1
    % update XData and YData   
    set(cart_1, 'Shape', cart_1_position{1, n});
    set(cart_2, 'Shape', cart_2_position{1, n});
    set(cart_3, 'Shape', cart_3_position{1, n});
    set(timePrint, 'String', sprintf('Time: %0.2f sec', timeAxis(n)));
    % get frame as an image
    frame = getframe(gcf);
    writeVideo(animation, frame);
end

close(animation);
end
