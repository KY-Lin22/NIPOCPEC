clear all
clc

%%
Data_NIP = load('Data_Newton_step_NIP');
Data_IP = load('Data_IP');
%%
continuationStepNum = Data_NIP.rec.Info{1}.continuationStepNum;
stepAxis = 0 : 1 : continuationStepNum;

figure(1)
semilogy(stepAxis, Data_NIP.rec.Info{1}.Log.VI_natural_residual, 'g*-',...
    stepAxis, Data_NIP.rec.Info{2}.Log.VI_natural_residual, 'bo-',...
    stepAxis, Data_NIP.rec.Info{3}.Log.VI_natural_residual, 'ks-',...
    stepAxis, Data_IP.rec.Info.Log.VI_natural_residual, 'rd-')
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
stairs(stepAxis, [Data_NIP.rec.Info{1}.Log.timeElapsed(2: end); Data_NIP.rec.Info{1}.Log.timeElapsed(end)], 'g-', 'LineWidth',1.1)
hold on
stairs(stepAxis, [Data_NIP.rec.Info{2}.Log.timeElapsed(2: end); Data_NIP.rec.Info{2}.Log.timeElapsed(end)], 'b-', 'LineWidth',1.1)
hold on
stairs(stepAxis, [Data_NIP.rec.Info{3}.Log.timeElapsed(2: end); Data_NIP.rec.Info{3}.Log.timeElapsed(end)], 'k-', 'LineWidth',1.1)
hold on
stairs(stepAxis, [Data_IP.rec.Info.Log.timeElapsed(2: end); Data_IP.rec.Info.Log.timeElapsed(end)], 'r--', 'LineWidth',1.1)
legend('NIP (1 corrector step)', 'NIP (2 corrector step)', 'NIP (3 corrector step)', 'IP', 'Location', 'northwest', 'FontSize', 11)
xlabel('Continuation step', 'FontSize', 11)    
ylabel('Time Elapsed [s]', 'FontSize', 11)
