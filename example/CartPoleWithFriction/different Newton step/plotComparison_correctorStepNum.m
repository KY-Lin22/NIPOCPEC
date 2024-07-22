clear all
clc

%%
Data_NIP = load('Data_Newton_step_NIP');
Data_IP = load('Data_IP');
%%
continuationStepNum = Data_NIP.rec.Info{1}.continuationStepNum;
stepAxis = 0 : 1 : continuationStepNum;

figure(1)
semilogy(stepAxis, Data_NIP.rec.Info{1}.Log.VI_nat_res, 'g*-',...
    stepAxis, Data_NIP.rec.Info{2}.Log.VI_nat_res, 'bo-',...
    stepAxis, Data_NIP.rec.Info{3}.Log.VI_nat_res, 'ks-',...
    stepAxis, Data_IP.rec.Info.Log.VI_nat_res, 'rd-')
grid on
legend('NIP (1 corrector step)', 'NIP (2 corrector step)', 'NIP (3 corrector step)', 'IP',...
    'Location','southwest', 'FontSize', 11)
xlabel('Continuation step', 'FontSize', 11)
ylabel('VI natural residual (max)', 'FontSize', 11)

figure(2)
semilogy(stepAxis, Data_NIP.rec.Info{1}.Log.KKT_error(:, end), 'g*-',...
    stepAxis, Data_NIP.rec.Info{2}.Log.KKT_error(:, end), 'bo-',...
    stepAxis, Data_NIP.rec.Info{3}.Log.KKT_error(:, end), 'ks-',...
    stepAxis, Data_IP.rec.Info.Log.KKT_error(:, end), 'rd-')

grid on
legend('NIP (1 corrector step)', 'NIP (2 corrector step)', 'NIP (3 corrector step)', 'IP', ...
    'Location', 'southwest', 'FontSize', 11)
xlabel('Continuation step', 'FontSize', 11)
ylabel('KKT total error', 'FontSize', 11)

%
figure(3)
stairs(stepAxis, [Data_NIP.rec.Info{1}.Log.time(2: end); Data_NIP.rec.Info{1}.Log.time(end)], 'g-', 'LineWidth',1.1)
hold on
stairs(stepAxis, [Data_NIP.rec.Info{2}.Log.time(2: end); Data_NIP.rec.Info{2}.Log.time(end)], 'b-', 'LineWidth',1.1)
hold on
stairs(stepAxis, [Data_NIP.rec.Info{3}.Log.time(2: end); Data_NIP.rec.Info{3}.Log.time(end)], 'k-', 'LineWidth',1.1)
hold on
stairs(stepAxis, [Data_IP.rec.Info.Log.time(2: end); Data_IP.rec.Info.Log.time(end)], 'r--', 'LineWidth',1.1)
grid on
set(gca, 'YScale', 'log')
legend('NIP (1 corrector step)', 'NIP (2 corrector step)', 'NIP (3 corrector step)', 'IP', 'Location', 'northwest', 'FontSize', 11)
xlabel('Continuation step', 'FontSize', 11)    
ylabel('Time Elapsed [s]', 'FontSize', 11)
