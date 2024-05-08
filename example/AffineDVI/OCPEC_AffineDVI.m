function OCPEC = OCPEC_AffineDVI()
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*
%%
% time parameter
timeHorizon = 1; % time horizon T
nStages = 100; % number of discretized stages
timeStep = timeHorizon ./ nStages; % discretization time step

% initial and reference state
x0 = [-0.5; -1]; % initial state
xRef = [0; 0]; % ref state

% variable and their bounds
xDim = 2;
uDim = 1;
lambdaDim = 1;
x = SX.sym('x', xDim, 1);
u = SX.sym('u', uDim, 1);
lambda = SX.sym('lambda', lambdaDim, 1);

xMax = [2; 2];
xMin = [-2; -2];
uMax = 2;
uMin = -2;

% cost function
xWeight = [20; 20];
uWeight = 1;
lambdaWeight = 1;
L_S = 0.5 * (x - xRef)'*diag(xWeight)*(x - xRef)...
    + 0.5 * u'*diag(uWeight)*u...
    + 0.5 * lambda'*diag(lambdaWeight)*lambda;
L_T = 0;

% DVI
f = [1, -3; -8, 10] * x + [4; 8] * u + [-3; -1] * lambda; % state equation f
F = [1, -3] * x + 3 * u + 5 * lambda; % VI function F
bl = -1;
bu = 1;
VISetType = 'box_constraint'; 

% inequality constraint G >= 0
G = [xMax - x;...
    x - xMin;...
    uMax - u;...
    u - uMin];
% equality constraint C = 0
C = SX(0,1);

%% create OCPEC instant
OCPEC = OCPEC_Formulation(...
    timeHorizon, nStages, timeStep,...
    x0,...
    x, u, lambda,...
    xMax, xMin, uMax, uMin,...
    L_T, L_S,...
    f, F, bl, bu, VISetType,...
    G, C);
end

