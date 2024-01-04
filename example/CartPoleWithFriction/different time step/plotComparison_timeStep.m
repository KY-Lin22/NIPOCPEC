clear all
clc

%%
Data_100ms_NIP = load('Data_T_3s_t_100ms_NIP');
Data_50ms_NIP = load('Data_T_3s_t_50ms_NIP');
Data_30ms_NIP = load('Data_T_3s_t_30ms_NIP');
Data_10ms_NIP = load('Data_T_3s_t_10ms_NIP');
Data_5ms_NIP = load('Data_T_3s_t_5ms_NIP');

Data_100ms_IP = load('Data_T_3s_t_100ms_IP');
Data_50ms_IP = load('Data_T_3s_t_50ms_IP');
Data_30ms_IP = load('Data_T_3s_t_30ms_IP');
Data_10ms_IP = load('Data_T_3s_t_10ms_IP');
Data_5ms_IP = load('Data_T_3s_t_5ms_IP');

Log_cost_100ms_NIP = Data_100ms_NIP.Info.Log.cost;
Log_cost_50ms_NIP = Data_50ms_NIP.Info.Log.cost;
Log_cost_30ms_NIP = Data_30ms_NIP.Info.Log.cost;
Log_cost_10ms_NIP = Data_10ms_NIP.Info.Log.cost;
Log_cost_5ms_NIP = Data_5ms_NIP.Info.Log.cost;

Log_cost_100ms_IP = Data_100ms_IP.Info.Log.cost;
Log_cost_50ms_IP = Data_50ms_IP.Info.Log.cost;
Log_cost_30ms_IP = Data_30ms_IP.Info.Log.cost;
Log_cost_10ms_IP = Data_10ms_IP.Info.Log.cost;
Log_cost_5ms_IP = Data_5ms_IP.Info.Log.cost;

continuationStepNum = Data_100ms_NIP.Info.continuationStepNum;

%% plot
index = 0 : 1 : continuationStepNum;
s_sequ = Data_100ms_NIP.Info.Log.param(:, 1);
figure(1)
semilogx(s_sequ, Log_cost_50ms_NIP, 'g*',...
    s_sequ, Log_cost_30ms_NIP, 'b*',...
    s_sequ, Log_cost_10ms_NIP, 'k*',...
    s_sequ, Log_cost_5ms_NIP, 'r*',...
    s_sequ, Log_cost_50ms_IP, 'go',...
    s_sequ, Log_cost_30ms_IP, 'bo',...
    s_sequ, Log_cost_10ms_IP, 'ko',...
    s_sequ, Log_cost_5ms_IP, 'ro')
grid on
label_50ms_NIP = '$\Delta t = 5 \cdot 10^{-2} (NIP)$';
label_30ms_NIP = '$\Delta t = 3 \cdot 10^{-2} (NIP)$';
label_10ms_NIP = '$\Delta t = 1 \cdot 10^{-2} (NIP)$';
label_5ms_NIP = '$\Delta t = 5 \cdot 10^{-3} (NIP)$';
label_50ms_IP = '$\Delta t = 5 \cdot 10^{-2} (IP)$';
label_30ms_IP = '$\Delta t = 3 \cdot 10^{-2} (IP)$';
label_10ms_IP = '$\Delta t = 1 \cdot 10^{-2} (IP)$';
label_5ms_IP = '$\Delta t = 5 \cdot 10^{-3} (IP)$';
legend(label_50ms_NIP, label_30ms_NIP, label_10ms_NIP, label_5ms_NIP,...
    label_50ms_IP, label_30ms_IP, label_10ms_IP, label_5ms_IP, ...
    'Interpreter','latex',...
    'Location','southwest')
xlabel('relaxation parameter s')
ylabel('cost')
ylim([600, 622])
