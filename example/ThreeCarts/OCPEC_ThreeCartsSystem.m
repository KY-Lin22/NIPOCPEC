function OCPEC = OCPEC_ThreeCartsSystem()
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*
%%
% time parameter
timeHorizon = 4; % time horizon T
nStages = 200; % number of discretized stages
timeStep = timeHorizon ./ nStages; % discretization time step

% physical parameter
m1 = 1; % cart mass
m2 = 1;
m3 = 1;
L1 = 2; % cart length
L2 = 2;
L3 = 2;
c = 2; % velocity damping ratio 

% initial state
x0 = [-3; 0; 3; 0; 0; 0]; 

% variable
xDim = 6; % position q1 q2 q3, velocity dq1 dq2 dq3
uDim = 1;
lambdaDim = 2;
x = SX.sym('x', xDim, 1);
u = SX.sym('u', uDim, 1);
lambda = SX.sym('lambda', lambdaDim, 1);

% cost function
xRef = [-5; 0; 5; 0; 0; 0]; % ref state
xWeight_S = [100; 100; 100; 0.1; 0.1; 0.1];
uWeight_S = 0.1;
lambdaWeight_S = 0.1;
L_S = 0.5 * (x - xRef)'*diag(xWeight_S)*(x - xRef)...
    + 0.5 * u'*diag(uWeight_S)*u ...
    + 0.5 * lambda'*diag(lambdaWeight_S)*lambda;
xWeight_T = [100; 10; 100; 10; 1; 10];
L_T = 0.5 * (x - xRef)'*diag(xWeight_T)*(x - xRef);
% L_T = 0;

% DVI
f = [x(4);...
    x(5);...
    x(6);...
    - (c/m1) * x(4)              - (1/m1) * lambda(1);...
    - (c/m2) * x(5) + (1/m2) * u + (1/m2) * lambda(1) - (1/m2) * lambda(2);...
    - (c/m3) * x(6)                                   + (1/m3) * lambda(2)];% state equation f
F = [x(2) - x(1) - 0.5 * (L1 + L2);...
    x(3) - x(2) - 0.5 * (L2 + L3)]; % VI function F
VISetType = 'nonnegative_orthant';
bl = zeros(lambdaDim, 1);
bu = inf * ones(lambdaDim, 1); 

% inequality constraint G >= 0
xMax = [8; 8; 8; 20; 20; 20];
xMin = [-8; -8; -8; -20; -20; -20];
uMax = 50;
uMin = -50;
G = [xMax - x;...
    x - xMin;...
    uMax - u;...
    u - uMin];
% equality constraint C = 0
C = SX(0,1);

%% create OCPEC instant
OCPEC = OCPEC_Formulation(...
    timeHorizon, nStages, timeStep,...
    x0, ...
    x, u, lambda,...
    xMax, xMin, uMax, uMin,...
    L_T, L_S,...
    f, F, bl, bu, VISetType,...
    G, C);
end