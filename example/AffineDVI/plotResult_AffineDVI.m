function plotResult_AffineDVI(OCPEC, NLP, z_Opt)
%UNTITLED Summary of this function goes here
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
plot(timeAxis, [OCPEC.x0(1), X_Opt(1, :)], 'r',...
     timeAxis, [OCPEC.x0(2), X_Opt(2, :)], 'g', 'LineWidth',1.2)
hold on
plot(timeAxis, -2*ones(1, numel(timeAxis)), 'k--',...
     timeAxis, 2*ones(1, numel(timeAxis)), 'k--',...
     'LineWidth', 0.5)
ylim([-2.3, 2.3])
legend('x1', 'x2')
xlabel('time [s]')
title('differential state')

subplot(3,1,2)
plot(timeAxis(2:end), U_Opt(1,:), 'LineWidth', 1.2)
hold on
plot(timeAxis, -2*ones(1, numel(timeAxis)), 'k--',...
     timeAxis, 2*ones(1, numel(timeAxis)), 'k--',...
     'LineWidth', 0.5)
ylim([-2.3, 2.3])
xlabel('time [s]')
title('control input')

subplot(3,1,3)
plot(timeAxis(2:end), LAMBDA_Opt(1, :), 'r',...
     timeAxis(2:end), F_Opt(1, :), 'b',...
     'LineWidth', 1.2)
hold on
plot(timeAxis, -1*ones(1, numel(timeAxis)), 'k--',...
     timeAxis, 1*ones(1, numel(timeAxis)), 'k--',...
     'LineWidth', 0.5)
ylim([1.3*min([LAMBDA_Opt(1, :), F_Opt(1, :)]),...
    1.3*max([LAMBDA_Opt(1, :), F_Opt(1, :)])])
legend('\lambda', 'F') 
xlabel('time [s]')
title('equilibrium constraint')
end