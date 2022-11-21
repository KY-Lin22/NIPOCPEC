clear all
clc

nStages = 300;
timeStep = 0.01;

option.mu = 1;
option.beta = 0.5;
option.T = 0.6;

InitPhase = [0; 0.5*pi; pi; 1.5*pi];
num_osci = size(InitPhase, 1);
InitRadius = option.mu * ones(num_osci, 1);
InitState = zeros(2 * num_osci, 1);
for i = 1 : num_osci
    InitState(1 + (i - 1) * 2 : i * 2, :) = [InitRadius(i)* cos(InitPhase(i)); InitRadius(i)* sin(InitPhase(i))];
end

[state, phase, radius] = holf_oscillator_couple(InitState, InitPhase, nStages, timeStep, option);

%%
osci_id = 1;
animateTraj_One_Oscillator(InitState(1 + (osci_id - 1) * 2 : osci_id * 2, :), InitPhase(1 + (osci_id - 1) * 2 : osci_id * 2, :),...
    state(1 + (osci_id - 1) * 2 : osci_id * 2, :), phase(1 + (osci_id - 1) * 2 : osci_id * 2, :), nStages, timeStep)

%%
timeAxis = 0 : timeStep : nStages * timeStep;

figure(1)
subplot(4,1,1)
plot(timeAxis, cos([InitPhase(1), phase(1,:)]), 'r', 'LineWidth',1.2)
legend('osci 1')
subplot(4,1,2)
plot(timeAxis, cos([InitPhase(2), phase(2,:)]), 'g', 'LineWidth',1.2)
legend('osci 2')
subplot(4,1,3)
plot(timeAxis, cos([InitPhase(3), phase(3,:)]), 'b', 'LineWidth',1.2)
legend('osci 3')
subplot(4,1,4)
plot(timeAxis, cos([InitPhase(4), phase(4,:)]), 'k', 'LineWidth',1.2)
legend('osci 4')
xlabel('time [s]')

