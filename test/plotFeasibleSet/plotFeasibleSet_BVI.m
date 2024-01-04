%% 
clear all
clc
s = 0.02;
bl = -1;
bu = 1;
plotFeasibleSet(s, bl, bu)


%%
function plotFeasibleSet(s, bl, bu)
% parameter
stepSize = 0.01;

% node
node_1_x = bl + 0.001; % bl upper node
node_1_y = s / (node_1_x - bl);

node_2_x = bl + 0.001; % bl lower node
node_2_y = - s / (bu - node_2_x);

node_3_x = bu - 0.001; % bu lower node
node_3_y = - s / (bu - node_3_x);

node_4_x = bu - 0.001; % bu upper node
node_4_y = s / (node_4_x - bl);

% region
curve_1_x = node_1_x : stepSize : node_4_x; % (lambda - bl) * eta = s
curve_1_y = zeros(1, length(curve_1_x));
for i = 1 : length(curve_1_x)
    curve_1_y(i) = s / (curve_1_x(i) - bl);
end
curve_2_x = node_3_x : -stepSize : node_2_x; % (bu - lambda) * eta = - s
curve_2_y = zeros(1, length(curve_2_x));
for i = 1 : length(curve_2_x)
    curve_2_y(i) = -s / (bu - curve_2_x(i));
end
reg_x = [curve_1_x, node_3_x, curve_2_x, node_1_x];
reg_y = [curve_1_y, node_3_y, curve_2_y, node_1_y];

% boundary
boundary_bl_y = node_2_y : stepSize : node_1_y;
boundary_bl_x = bl * ones(1, length(boundary_bl_y));
boundary_bu_y = node_3_y : stepSize : node_4_y;
boundary_bu_x = bu * ones(1, length(boundary_bu_y));

%% plot
figure(1)
% relaxed region
patch(reg_x, reg_y, 'blue', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
% boundary
plot(boundary_bl_x, boundary_bl_y, 'k', 'LineWidth', 2)
hold on
plot(boundary_bu_x, boundary_bu_y, 'k', 'LineWidth', 2)
hold on
plot(curve_1_x, curve_1_y, 'k', 'LineWidth', 2)
hold on
plot(curve_2_x, curve_2_y, 'k', 'LineWidth', 2)
hold on
axis([bl - 0.4, bu + 0.4, node_2_y - 1, node_4_y + 1])
xline(0, 'LineWidth', 1); 
yline(0, 'LineWidth', 1);

end