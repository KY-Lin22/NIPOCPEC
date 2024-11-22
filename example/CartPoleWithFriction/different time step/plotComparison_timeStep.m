clear all
clc

%%
Data_NIP = load('Data_time_step_NIP');
Data_IP = load('Data_time_step_IP');
Data_IP_single = load('Data_time_step_IP_single');

s_sequence = Data_NIP.rec.Info{1}.Log.param(:, 1);
timeHorizon = Data_NIP.rec.timeHorizon;
nStage_sequence = Data_NIP.rec.nStage_sequence;
timeStep = [...
    timeHorizon/nStage_sequence{1};...
    timeHorizon/nStage_sequence{2};...
    timeHorizon/nStage_sequence{3}];
axis_limit_bais = 200;

%% cold start
figure(1)
subplot(3, 1, 1)
IP_single_1_cost = Data_IP_single.rec.Info{1}.Log.cost;
y_min_cold_1 = min(IP_single_1_cost) - axis_limit_bais;
y_max_cold_1 = max(IP_single_1_cost) + axis_limit_bais;
semilogx(s_sequence, IP_single_1_cost, 'bx');
hold on
semilogx([timeStep(1), timeStep(1)], [y_min_cold_1, y_max_cold_1], 'k--', 'LineWidth', 1.2);
grid on
title('$\Delta t = 5 \cdot 10^{-2}$', 'Interpreter','latex', 'FontSize', 11)
ylabel('Cost', 'FontSize', 11)
ylim([y_min_cold_1, y_max_cold_1])

subplot(3, 1, 2)
IP_single_2_cost = Data_IP_single.rec.Info{2}.Log.cost;
y_min_cold_2 = min(IP_single_2_cost) - axis_limit_bais;
y_max_cold_2 = max(IP_single_2_cost) + axis_limit_bais;

semilogx(s_sequence, IP_single_2_cost, 'bx');
hold on
semilogx([timeStep(2), timeStep(2)], [y_min_cold_2, y_max_cold_2], 'k--', 'LineWidth', 1.2);
grid on
title('$\Delta t = 1 \cdot 10^{-2}$', 'Interpreter','latex', 'FontSize', 11)
ylabel('Cost', 'FontSize', 11)
ylim([y_min_cold_2, y_max_cold_2])

subplot(3, 1, 3)
IP_single_3_cost = Data_IP_single.rec.Info{3}.Log.cost;
y_min_cold_3 = min(IP_single_3_cost) - axis_limit_bais;
y_max_cold_3 = max(IP_single_3_cost) + axis_limit_bais;

semilogx(s_sequence, IP_single_3_cost, 'bx');
hold on
semilogx([timeStep(3), timeStep(3)], [y_min_cold_3, y_max_cold_3], 'k--', 'LineWidth', 1.2);
grid on
title('$\Delta t = 5 \cdot 10^{-3}$', 'Interpreter','latex', 'FontSize', 11)
ylabel('Cost', 'FontSize', 11)
ylim([y_min_cold_3, y_max_cold_3])
xlabel('Relaxation parameter', 'FontSize', 11)

%% continuation method
figure(2)
subplot(3, 1, 1)
NIP_1_cost = Data_NIP.rec.Info{1}.Log.cost;
IP_1_cost = Data_IP.rec.Info{1}.Log.cost;
y_min_cont_1 = min([min(NIP_1_cost), min(IP_1_cost)]) - 8;
y_max_cont_1 = max([max(NIP_1_cost), max(IP_1_cost)]) + 3;

semilogx(s_sequence, NIP_1_cost, 'g*');
hold on
semilogx(s_sequence, IP_1_cost, 'ro');
hold on
semilogx([timeStep(1), timeStep(1)], [y_min_cont_1, y_max_cont_1], 'k--', 'LineWidth', 1.2);
grid on
title('$\Delta t = 5 \cdot 10^{-2}$', 'Interpreter','latex', 'FontSize', 11)
ylabel('Cost', 'FontSize', 11)
ylim([y_min_cont_1, y_max_cont_1])

subplot(3, 1, 2)
NIP_2_cost = Data_NIP.rec.Info{2}.Log.cost;
IP_2_cost = Data_IP.rec.Info{2}.Log.cost;
y_min_cont_2 = min([min(NIP_2_cost), min(IP_2_cost)]) - 1;
y_max_cont_2 = max([max(NIP_2_cost), max(IP_2_cost)]) + 1;

semilogx(s_sequence, NIP_2_cost, 'g*');
hold on
semilogx(s_sequence, IP_2_cost, 'ro');
hold on
semilogx([timeStep(2), timeStep(2)], [y_min_cont_2, y_max_cont_2], 'k--', 'LineWidth', 1.2);
grid on
title('$\Delta t = 1 \cdot 10^{-2}$', 'Interpreter','latex', 'FontSize', 11)
ylabel('Cost', 'FontSize', 11)
ylim([y_min_cont_2, y_max_cont_2])

subplot(3, 1, 3)
NIP_3_cost = Data_NIP.rec.Info{3}.Log.cost;
IP_3_cost = Data_IP.rec.Info{3}.Log.cost;
y_min_cont_3 = min([min(NIP_3_cost), min(IP_3_cost)]) - 1;
y_max_cont_3 = max([max(NIP_3_cost), max(IP_3_cost)]) + 1;

semilogx(s_sequence, NIP_3_cost, 'g*');
hold on
semilogx(s_sequence, IP_3_cost, 'ro');
hold on
semilogx([timeStep(3), timeStep(3)], [y_min_cont_3, y_max_cont_3], 'k--', 'LineWidth', 1.2);
grid on
title('$\Delta t = 5 \cdot 10^{-3}$', 'Interpreter','latex', 'FontSize', 11)
xlabel('Relaxation parameter', 'FontSize', 11)
ylabel('Cost', 'FontSize', 11)
ylim([y_min_cont_3, y_max_cont_3])