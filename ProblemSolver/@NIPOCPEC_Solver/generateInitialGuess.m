function generateInitialGuess(solver)
%generateInitialGuess
%   generate Initial Guess

OCPEC = solver.OCPEC;
nStages = OCPEC.nStages;
Dim = OCPEC.Dim;
plant = OCPEC.plant;
% init
Iterate = struct('tau', zeros(Dim.tau, nStages), 'x', zeros(Dim.x, nStages),...
    'p', zeros(Dim.p, nStages), 'w', zeros(Dim.w, nStages),...
    'sigma', zeros(Dim.sigma, nStages), 'eta', zeros(Dim.eta, nStages),...
    'lambda', zeros(Dim.lambda, nStages), 'gamma', zeros(Dim.gamma, nStages));

disp('Generating Initial Guess...')

%% generate initial guess for primal variable
% tau
tau_Init = randn(Dim.tau, 1);
tau_End = OCPEC.TerminalCost.tauRef;
Iterate.tau = TrajectoryInterpolation(tau_Init, tau_End, nStages);
% x
x_Init = OCPEC.InitState;
x_End = OCPEC.TerminalCost.xRef;
Iterate.x = TrajectoryInterpolation(x_Init, x_End, nStages);
% p
Iterate.p = zeros(Dim.p, nStages);
for i = 1 : Dim.p
    if (plant.l(i) == 0) && (plant.u(i) == Inf)
        % nonlinear complementary problem
        Iterate.p(i, :) = abs(randn(1, nStages)); % p > = 0
    else
        % box constraint variation inequality
        Iterate.p(i, :) = repmat(1/2*(plant.l(i) + plant.u(i)), 1, nStages); % l < = p < = u
    end
end
% w
if (strcmp(OCPEC.VI_mode, 'SmoothingEquation')) || (strcmp(OCPEC.VI_mode, 'Reg_Scholtes'))
    Iterate.w = zeros(Dim.w, nStages);
    % w = K
    for n = 1 : nStages
        Iterate.w(:, n) = plant.computeVIFunc(Iterate.tau(:, n), Iterate.x(:, n), Iterate.p(:, n)); 
    end   
elseif strcmp(OCPEC.VI_mode, 'Reg_NCPs')
    % w >= 0
    Iterate.w = ones(Dim.w, nStages); 
end

%% generate initial guess for dual variable
Iterate.sigma  = ones(Dim.sigma, nStages); % sigma >=0
Iterate.eta    = randn(Dim.eta, nStages);
Iterate.lambda = randn(Dim.lambda, nStages);
if (strcmp(OCPEC.VI_mode, 'Reg_NCPs')) || (strcmp(OCPEC.VI_mode, 'Reg_Scholtes'))
    Iterate.gamma  = ones(Dim.gamma, nStages); % gamma >=0
elseif (strcmp(OCPEC.VI_mode, 'SmoothingEquation'))
    Iterate.gamma  = randn(Dim.gamma, nStages);
end

%% save initial guess
save('Gen_InitialGuess.mat', 'Iterate');
disp('Done!')

end

