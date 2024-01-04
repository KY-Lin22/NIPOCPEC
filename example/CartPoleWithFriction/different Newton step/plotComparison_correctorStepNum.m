clear all
clc

%%
Data_corrector_1step_NIP = load('Data_corrector_1step_NIP');
Data_corrector_2step_NIP = load('Data_corrector_2step_NIP');
Data_corrector_3step_NIP = load('Data_corrector_3step_NIP');
Data_IP = load('Data_IP');

Log_KKT_error_primal_1step_NIP = Data_corrector_1step_NIP.Info.Log.KKT_error(:, 1);
Log_KKT_error_primal_2step_NIP = Data_corrector_2step_NIP.Info.Log.KKT_error(:, 1);
Log_KKT_error_primal_3step_NIP = Data_corrector_3step_NIP.Info.Log.KKT_error(:, 1);
Log_KKT_error_primal_IP = Data_IP.Info.Log.KKT_error(:, 1);

Log_KKT_error_dual_1step_NIP = Data_corrector_1step_NIP.Info.Log.KKT_error(:, 3);
Log_KKT_error_dual_2step_NIP = Data_corrector_2step_NIP.Info.Log.KKT_error(:, 3);
Log_KKT_error_dual_3step_NIP = Data_corrector_3step_NIP.Info.Log.KKT_error(:, 3);
Log_KKT_error_dual_IP = Data_IP.Info.Log.KKT_error(:, 2);

Log_KKT_error_total_1step_NIP = Data_corrector_1step_NIP.Info.Log.KKT_error(:, 6);
Log_KKT_error_total_2step_NIP = Data_corrector_2step_NIP.Info.Log.KKT_error(:, 6);
Log_KKT_error_total_3step_NIP = Data_corrector_3step_NIP.Info.Log.KKT_error(:, 6);
Log_KKT_error_total_IP = Data_IP.Info.Log.KKT_error(:, 3);


Log_VI_natRes_1step_NIP = Data_corrector_1step_NIP.Info.Log.VI_natural_residual;
Log_VI_natRes_2step_NIP = Data_corrector_2step_NIP.Info.Log.VI_natural_residual;
Log_VI_natRes_3step_NIP = Data_corrector_3step_NIP.Info.Log.VI_natural_residual;
Log_VI_natRes_IP = Data_IP.Info.Log.VI_natural_residual;

Log_timeElapsed_1step_NIP = Data_corrector_1step_NIP.Info.Log.timeElapsed;
Log_timeElapsed_2step_NIP = Data_corrector_2step_NIP.Info.Log.timeElapsed;
Log_timeElapsed_3step_NIP = Data_corrector_3step_NIP.Info.Log.timeElapsed;
Log_timeElapsed_IP = Data_IP.Info.Log.timeElapsed;

continuationStepNum = Data_corrector_1step_NIP.Info.continuationStepNum;

%%
stepAxis = 0 : 1 : continuationStepNum;
% figure(1)
% semilogy(index, Log_KKT_error_primal_1step_NIP, 'g*-',...
%     index, Log_KKT_error_primal_2step_NIP, 'bo-',...
%     index, Log_KKT_error_primal_3step_NIP, 'ks-',...
%     index, Log_KKT_error_primal_IP, 'rd-')
% grid on
% legend('NIP (1 step)', 'NIP (2 step)', 'NIP (3 step)', 'IP')
% xlabel('continuation step')
% ylabel('KKT primal error')

% figure(2)
% semilogy(index, Log_KKT_error_dual_1step_NIP, 'g*-',...
%     index, Log_KKT_error_dual_2step_NIP, 'bo-',...
%     index, Log_KKT_error_dual_3step_NIP, 'ks-',...
%     index, Log_KKT_error_dual_IP, 'rd-')
% grid on
% legend('NIP (1 step)', 'NIP (2 step)', 'NIP (3 step)', 'IP')
% xlabel('continuation step')
% ylabel('KKT dual error')

figure(3)
semilogy(stepAxis, Log_KKT_error_total_1step_NIP, 'g*-',...
    stepAxis, Log_KKT_error_total_2step_NIP, 'bo-',...
    stepAxis, Log_KKT_error_total_3step_NIP, 'ks-',...
    stepAxis, Log_KKT_error_total_IP, 'rd-')
grid on
legend('NIP (1 corrector step)', 'NIP (2 corrector step)', 'NIP (3 corrector step)', 'IP' , 'Location','southwest')
xlabel('continuation step')
ylabel('KKT total error')

figure(4)
semilogy(stepAxis, Log_VI_natRes_1step_NIP, 'g*-',...
    stepAxis, Log_VI_natRes_2step_NIP, 'bo-',...
    stepAxis, Log_VI_natRes_3step_NIP, 'ks-',...
    stepAxis, Log_VI_natRes_IP, 'rd-')
grid on
legend('NIP (1 corrector step)', 'NIP (2 corrector step)', 'NIP (3 corrector step)', 'IP' , 'Location','southwest')
xlabel('continuation step')
ylabel('VI natural residual (max)')

figure(5)
stairs(stepAxis, [Log_timeElapsed_1step_NIP(2: end); Log_timeElapsed_1step_NIP(end)], 'g-', 'LineWidth',1.1)
hold on
stairs(stepAxis, [Log_timeElapsed_2step_NIP(2: end); Log_timeElapsed_2step_NIP(end)], 'b-', 'LineWidth',1.1)
hold on
stairs(stepAxis, [Log_timeElapsed_3step_NIP(2: end); Log_timeElapsed_3step_NIP(end)], 'k-', 'LineWidth',1.1)
hold on
stairs(stepAxis, [Log_timeElapsed_IP(2: end); Log_timeElapsed_IP(end)], 'r--', 'LineWidth',1.1)
legend('NIP (1 corrector step)', 'NIP (2 corrector step)', 'NIP (3 corrector step)', 'IP', 'Location','northwest')
xlabel('continuation step')    
ylabel('time Elapsed [s]')
