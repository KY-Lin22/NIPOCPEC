function OCPEC = OCPEC_Acrobot()
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*

% time parameter
timeHorizon = 5; % time horizon T
nStages = 100; % number of discretized stages
timeStep = timeHorizon ./ nStages; % discretization time step

% physical parameter
m1 = 1; % mass
m2 = 0.1;
L1 = 1; % link length
L2 = 1; 
r1 = 0.5; % link center
r2 = 0.5;
I1 = 0.33; % inertia
I2 = 0.33;
JF1 = 0.1; % joint friction
JF2 = 0.1;
q2_min = -1/2*pi; % joint 2 max and min angle
q2_max = 1/2*pi;
gravity_param = 9.8; % gravity param
% initial state
x0 = [0/180*pi; 0/180*pi; 0; 0];

% variable
xDim = 4; % position q1 q2, velocity dq1 dq2
uDim = 1;
lambdaDim = 2;
x = SX.sym('x', xDim, 1);
u = SX.sym('u', uDim, 1);
lambda = SX.sym('lambda', lambdaDim, 1);

% cost function
xRef = [180/180*pi; 0/180*pi; 0; 0];
xWeight_S = [5; 5; 0.01; 0.01];
uWeight_S = 0.01;
lambdaWeight_S = 0.0001;
L_S = 0.5 * (x - xRef)'*diag(xWeight_S)*(x - xRef)...
    + 0.5 * u'*diag(uWeight_S)*u ...
    + 0.5 * lambda'*diag(lambdaWeight_S)*lambda;
xWeight_T = [5; 5; 0.1; 0.1];
L_T = 0.5 * (x - xRef)'*diag(xWeight_T)*(x - xRef);

% DVI
M_Matrix = [I1 + I2 + m2*L1^2 + 2*m2*L1*r2*cos(x(2)), I2 + m2*L1*r2*cos(x(2));...
            I2 + m2*L1*r2*cos(x(2)),                  I2]; % mass
C_Matrix = [-2*m2*L1*r2*sin(x(2))*x(4),    -m2*L1*r2*sin(x(2))*x(4);...
            m2*L1*r2*sin(x(2))*x(3),        0]; % nonlinear terms
G_Matrix = [-m1*gravity_param*r1*sin(x(1)) - m2*gravity_param*(L1*sin(x(1)) + r2*sin(x(1)+x(2)));...
            -m2*gravity_param*r2*sin(x(1)+x(2))]; % gravity
Bu_Matrix = [0;...
            u(1)];% control
P_Matrix = [0;...
            lambda(1) - lambda(2)]; % external force
H_Matrix = G_Matrix + Bu_Matrix + P_Matrix - C_Matrix * [x(3); x(4)] - [JF1 * x(3); JF2 * x(4)];  

f = [x(3);...
     x(4);...
     inv(M_Matrix)*H_Matrix];

F = [x(2) - q2_min;...
     q2_max - x(2)];
VISetType = 'nonnegative_orthant';
bl = zeros(lambdaDim, 1);
bu = inf * ones(lambdaDim, 1); 
% inequality constraint G >= 0
xMax = [2*pi; q2_max; 30; 30];
xMin = [-2*pi; q2_min; -30; -30];
uMax = 30;
uMin = -30;
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