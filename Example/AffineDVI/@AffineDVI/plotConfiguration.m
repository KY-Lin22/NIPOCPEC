function plotConfiguration(plant, InitState, RefState)
%plotConfiguration
%   Detailed explanation goes here

DynVarLimit = plant.DynVarLimit;
baseLine_horizontal_X = [DynVarLimit.x_Min(1); DynVarLimit.x_Max(1)];
baseLine_horizontal_Y = [0; 0];  
baseLine_vertical_X = [0; 0];
baseLine_vertical_Y = [DynVarLimit.x_Min(2); DynVarLimit.x_Max(2)]; 

figure(2000);
plot(baseLine_horizontal_X, baseLine_horizontal_Y, '.-k', 'MarkerSize', 0.5, 'LineWidth', 0.5);
hold on
plot(baseLine_vertical_X, baseLine_vertical_Y, '.-k', 'MarkerSize', 0.5, 'LineWidth', 0.5);
hold on
plot(InitState(1), InitState(2), '*', 'Color', [1 0 0], 'MarkerSize', 6, 'LineWidth', 2);
hold on
plot(RefState(1), RefState(2),  'o', 'Color', [1 0 0], 'MarkerSize', 6, 'LineWidth', 2);
legend('Init', 'Ref')
axis equal
axisLimit_X = baseLine_horizontal_X;
axisLimit_Y = baseLine_vertical_Y; 
axis([axisLimit_X; axisLimit_Y]); 

end

