clear all
clc
%% original feasible region
x_Ori = [0, 0, 1];
y_Ori = [1, 0, 0];
plot(x_Ori, y_Ori, 'k', 'LineWidth', 5)
axis([-0.2 1 -0.2 1])
xlabel('g_i')
ylabel('\eta_i')

%% penalty method
% original
x_Ori = [0, 0, 1];
y_Ori = [1, 0, 0];
% penalty
x_Pen = [0, 1, 1, 0];
y_Pen = [0, 0, 1, 1];

figure(1)
patch(x_Pen,y_Pen,'green')
hold on
plot(x_Ori, y_Ori, 'k', 'LineWidth', 5)
axis([-0.2 1 -0.2 1])
xlabel('g_i')
ylabel('\eta_i')

%% smoothing method
% original
x_Ori = [0, 0, 1];
y_Ori = [1, 0, 0];
% smoothing
s = 0.01;
x_Init = 0.01;
stepsize = 0.01;
x_End = 1;
x_intepolation = x_Init: stepsize : x_End;
y_intepolation = s ./ x_intepolation;

figure(2)
plot(x_intepolation, y_intepolation, 'g', 'LineWidth', 3)
hold on
plot(x_Ori, y_Ori, 'k', 'LineWidth', 5)
axis([0 1 0 1])
xlabel('g_i')
ylabel('\eta_i')

%% relaxation method
% original
x_Ori = [0, 0, 1];
y_Ori = [1, 0, 0];
% relaxation
s = 0.01;
s_1 = 0.05;
s_2 = 0.1;
x_Init = 0.01;
stepsize = 0.01;
x_End = 1;
x_intepolation = x_Init: stepsize : x_End;
y_intepolation = s ./ x_intepolation;
x_relaxation = [x_intepolation, x_intepolation(1), x_intepolation(1)];
y_relaxation = [y_intepolation, y_intepolation(end), y_intepolation(1)];

y_intepolation_1 = s_1 ./ x_intepolation;
y_intepolation_2 = s_2 ./ x_intepolation;

figure(3)
patch(x_relaxation,y_relaxation,'green')
hold on
plot(x_Ori, y_Ori, 'k', 'LineWidth', 5)
hold on
plot(x_intepolation, y_intepolation_1)
hold on
plot(x_intepolation, y_intepolation_2)
axis([-0.2 1 -0.2 1])
xlabel('g_i')
ylabel('\eta_i')