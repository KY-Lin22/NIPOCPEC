function plotResult_ThreeCartsSystem(OCPEC, NLP, z_Opt)
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
     timeAxis, [OCPEC.x0(2), X_Opt(2, :)], 'g',...
     timeAxis, [OCPEC.x0(3), X_Opt(3, :)], 'b', 'LineWidth',1.2)
legend('cart 1', 'cart 2', 'cart 3')
xlabel('time [s]')
ylabel('position [m]')
title('cart position')

subplot(3,1,2)
plot(timeAxis, [OCPEC.x0(4), X_Opt(4, :)], 'r',...
     timeAxis, [OCPEC.x0(5), X_Opt(5, :)], 'g',...
     timeAxis, [OCPEC.x0(6), X_Opt(6, :)], 'b', 'LineWidth',1.2)
legend('cart 1', 'cart 2', 'cart 3')
xlabel('time [s]')
ylabel('velocity [m/s]')
title('cart velocity')

subplot(3,1,3)
plot(timeAxis(2:end), U_Opt(1,:), 'LineWidth', 1.2)
xlabel('time [s]')
title('actuated force')

figure(112)
subplot(2,1,1)
plot(timeAxis(2:end), F_Opt(1, :), 'k', ...
    timeAxis(2:end), F_Opt(2, :), 'b','LineWidth', 1.2)
legend('cart 1-2', 'cart 2-3 ') 
xlabel('time [s]')
title('cart distance')

subplot(2,1,2)
plot(timeAxis(2:end), LAMBDA_Opt(1, :), 'k',...
    timeAxis(2:end), LAMBDA_Opt(2, :), 'b', 'LineWidth', 1.2)
legend('cart 1-2', 'cart 2-3 ') 
xlabel('time [s]')
title('contact force')

end