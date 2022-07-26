function showInfo(OCPEC)
%showInfo
%   Detailed explanation goes here

%% show OCPEC problem information
%
Dim = OCPEC.Dim;
nStages = OCPEC.nStages;
timeStep = OCPEC.timeStep;
disp('*----------------- Perturbed OCPEC Problem Information ----------------------*')
disp(['Perturbed Reformulation of Variational Inequality: ', OCPEC.VI_mode])
% time
disp(['timeHorizon: ...', num2str(nStages * timeStep), ' s'])
disp(['- Stage: .......', num2str(nStages),])
disp(['- timeStep: ....', num2str(timeStep),' s '])
% number of variable
disp(['Number of Total Variables: .................', num2str(Dim.Y),   ' * ', num2str(nStages), ' = ', num2str(Dim.Y * nStages)])
disp(['Number of Primal Variable: .................', num2str(Dim.Z),   ' * ', num2str(nStages), ' = ', num2str(Dim.Z * nStages)])
disp(['- tau(control input): ......................', num2str(Dim.tau), ' * ', num2str(nStages), ' = ', num2str(Dim.tau * nStages), '; '])
disp(['- x(system state): .........................', num2str(Dim.x),   ' * ', num2str(nStages), ' = ', num2str(Dim.x * nStages), '; '])
disp(['- p(equilibrium state): ....................', num2str(Dim.p),   ' * ', num2str(nStages), ' = ', num2str(Dim.p * nStages), '; '])
disp(['- w(auxiliary variable): ...................', num2str(Dim.w),   ' * ', num2str(nStages), ' = ', num2str(Dim.w * nStages), '; '])

disp(['Number of Dual Variable: ...................',  num2str(Dim.LAMBDA), ' * ', num2str(nStages), ' = ', num2str(Dim.LAMBDA * nStages)])
disp(['- sigma(inequality): .......................', num2str(Dim.sigma),   ' * ', num2str(nStages), ' = ', num2str(Dim.sigma * nStages), '; '])
disp(['- eta(equality): ...........................', num2str(Dim.eta),     ' * ', num2str(nStages), ' = ', num2str(Dim.eta * nStages), '; '])
disp(['- lambda(state equation): ..................', num2str(Dim.lambda),  ' * ', num2str(nStages), ' = ', num2str(Dim.lambda * nStages), '; '])
disp(['- gamma(perturbed equilibrium dynamics): ...', num2str(Dim.gamma),   ' * ', num2str(nStages), ' = ', num2str(Dim.gamma * nStages)])
disp([])
disp('*----------------------------------------------------------------------------*')

end

