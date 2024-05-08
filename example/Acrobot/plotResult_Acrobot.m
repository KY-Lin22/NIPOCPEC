function plotResult_Acrobot(OCPEC, NLP, z_Opt)
%UNTITLED21 Summary of this function goes here
%   Detailed explanation goes here
Z_Opt = reshape(z_Opt, NLP.Dim.z_Node(end), OCPEC.nStages);

X_Opt = Z_Opt(1 : NLP.Dim.z_Node(1), :);
U_Opt = Z_Opt(NLP.Dim.z_Node(1) + 1 : NLP.Dim.z_Node(2), :);
LAMBDA_Opt = Z_Opt(NLP.Dim.z_Node(2) + 1 : NLP.Dim.z_Node(3), :);

timeAxis = 0 : OCPEC.timeStep : OCPEC.nStages * OCPEC.timeStep;

F_FuncObj_map = OCPEC.FuncObj.F.map(OCPEC.nStages);
F_Opt = full(F_FuncObj_map(X_Opt, U_Opt, LAMBDA_Opt));

figure(111)
subplot(3,1,1)
plot(timeAxis, [OCPEC.x0(1), X_Opt(1, :)]*360/(2*pi), 'r',...
     timeAxis, [OCPEC.x0(2), X_Opt(2, :)]*360/(2*pi), 'g', 'LineWidth',1.2)
legend('link 1', 'link 2')
xlabel('time [s]')
ylabel('angle [deg]')
title('link angle')

subplot(3,1,2)
plot(timeAxis, [OCPEC.x0(3), X_Opt(3, :)], 'r',...
     timeAxis, [OCPEC.x0(4), X_Opt(4, :)], 'g', 'LineWidth',1.2)
legend('link 1', 'link 2')
xlabel('time [s]')
ylabel('velocity [rad/s]')
title('link velocity')

subplot(3,1,3)
plot(timeAxis(2:end), U_Opt(1,:), 'LineWidth', 1.2)
xlabel('time [s]')
title('actuated force')

figure(112)
subplot(2,1,1)
plot(timeAxis(2:end), F_Opt(1, :), 'k', ...
    timeAxis(2:end), F_Opt(2, :), 'b','LineWidth', 1.2)
legend('q2 - qmin', 'qmax - q2') 
xlabel('time [s]')
title('joint limit distance')

subplot(2,1,2)
plot(timeAxis(2:end), LAMBDA_Opt(1, :), 'k',...
    timeAxis(2:end), LAMBDA_Opt(2, :), 'b', 'LineWidth', 1.2)
legend('qmin force', 'qmax force') 
xlabel('time [s]')
title('contact force')

end