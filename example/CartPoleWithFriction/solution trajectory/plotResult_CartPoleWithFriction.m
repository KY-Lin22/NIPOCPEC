function plotResult_CartPoleWithFriction(OCPEC, NLP, z_Opt)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
Z_Opt = reshape(z_Opt, NLP.Dim.z_Node(end), OCPEC.nStages);

X_Opt = Z_Opt(1 : NLP.Dim.z_Node(1), :);
U_Opt = Z_Opt(NLP.Dim.z_Node(1) + 1 : NLP.Dim.z_Node(2), :);
LAMBDA_Opt = Z_Opt(NLP.Dim.z_Node(2) + 1 : NLP.Dim.z_Node(3), :);

timeAxis = 0 : OCPEC.timeStep : OCPEC.nStages * OCPEC.timeStep;

F_FuncObj_map = OCPEC.FuncObj.F.map(OCPEC.nStages);
F_Opt = full(F_FuncObj_map(X_Opt, U_Opt, LAMBDA_Opt));

figure(111)
subplot(4,1,1)
plot(timeAxis, [OCPEC.x0(1), X_Opt(1, :)], 'r',...
     timeAxis, [OCPEC.x0(2), X_Opt(2, :)], 'g', 'LineWidth',1.2)
legend('$x_c$', '$\theta_p$',  'Interpreter','latex', 'FontSize', 11, 'Location', 'northeast')
legend('boxoff')
%xlabel('time [s]')
ylim([-2, 4])
ylabel('Position', 'FontSize', 11) 

subplot(4,1,2)
plot(timeAxis, [OCPEC.x0(3), X_Opt(3, :)], 'r',...
     timeAxis, [OCPEC.x0(4), X_Opt(4, :)], 'g', 'LineWidth',1.2)
legend('$\dot{x}_c$', '$\dot{\theta}_p$', 'Interpreter','latex', 'FontSize', 11, 'Location', 'northeast')
legend('boxoff')
%xlabel('time [s]')
ylim([-6, 8])
ylabel('Velocity', 'FontSize', 11)  

subplot(4,1,3)
plot(timeAxis(2:end), U_Opt(1,:), 'LineWidth', 1.2)
legend('$u$', 'Interpreter', 'latex', 'FontSize', 11)
legend('boxoff')
%xlabel('time [s]')
ylim([-35, 35])
ylabel('Control', 'FontSize', 11)

subplot(4,1,4)
plot(timeAxis(2:end), LAMBDA_Opt(1, :), 'r',...
     timeAxis(2:end), F_Opt(1, :), 'b',...
     'LineWidth', 1.2)
legend('$\lambda$', '$F$', 'Interpreter','latex', 'FontSize', 11, 'Location', 'southeast') 
legend('boxoff')
xlabel('Time [s]', 'FontSize', 11)
ylim([-6, 4])
ylabel('Friction', 'FontSize', 11)
end

