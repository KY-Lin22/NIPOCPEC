function z_Init = preprocessInitialGuess(self, z_Init)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*
% load parameter
xDim = self.OCPEC.Dim.x;
uDim = self.OCPEC.Dim.u;
lambdaDim = self.OCPEC.Dim.lambda;
nStages = self.OCPEC.nStages;
xMin = self.OCPEC.xMin;
xMax = self.OCPEC.xMax;
uMin = self.OCPEC.uMin;
uMax = self.OCPEC.uMax;
lambdaMin = self.OCPEC.bl;
lambdaMax = self.OCPEC.bu;

kappa_bound = self.Option.Init.kappa_bound;
kappa_interval = self.Option.Init.kappa_interval;

%% proprocess
% extract x, u, lambda, eta
Z_Init = reshape(z_Init, self.NLP.Dim.z_Node(end), nStages);

X_Init      = Z_Init(                         1 : self.NLP.Dim.z_Node(1), :);
U_Init      = Z_Init(self.NLP.Dim.z_Node(1) + 1 : self.NLP.Dim.z_Node(2), :);
LAMBDA_Init = Z_Init(self.NLP.Dim.z_Node(2) + 1 : self.NLP.Dim.z_Node(3), :);

% x_Init
X_lb = xMin + min([kappa_bound * max([ones(xDim, 1), abs(xMin)], [], 2), kappa_interval * (xMax - xMin)], [], 2);
X_ub = xMax - min([kappa_bound * max([ones(xDim, 1), abs(xMax)], [], 2), kappa_interval * (xMax - xMin)], [], 2);
x_lb = repmat(X_lb, nStages, 1);
x_ub = repmat(X_ub, nStages, 1);
x_Init = reshape(X_Init, [], 1);
x_Init = min([max([x_lb, x_Init], [], 2), x_ub], [], 2); % project into [x_lb, x_ub]
X_Init = reshape(x_Init, xDim, nStages);

% u_Init
U_lb = uMin + min([kappa_bound * max([ones(uDim, 1), abs(uMin)], [], 2), kappa_interval * (uMax - uMin)], [], 2);
U_ub = uMax - min([kappa_bound * max([ones(uDim, 1), abs(uMax)], [], 2), kappa_interval * (uMax - uMin)], [], 2);
u_lb = repmat(U_lb, nStages, 1);
u_ub = repmat(U_ub, nStages, 1);
u_Init = reshape(U_Init, [], 1);
u_Init = min([max([u_lb, u_Init], [], 2), u_ub], [], 2); % project into [u_lb, u_ub]
U_Init = reshape(u_Init, uDim, nStages);

% lambda_Init
LAMBDA_lb = lambdaMin + min([kappa_bound * max([ones(lambdaDim, 1), abs(lambdaMin)], [], 2), kappa_interval * (lambdaMax - lambdaMin)], [], 2);
LAMBDA_ub = lambdaMax - min([kappa_bound * max([ones(lambdaDim, 1), abs(lambdaMax)], [], 2), kappa_interval * (lambdaMax - lambdaMin)], [], 2);
lambda_lb = repmat(LAMBDA_lb, nStages, 1);
lambda_ub = repmat(LAMBDA_ub, nStages, 1);
lambda_Init = reshape(LAMBDA_Init, [], 1);
lambda_Init = min([max([lambda_lb, lambda_Init], [], 2), lambda_ub], [], 2); % project into [lambda_lb, lambda_ub]
LAMBDA_Init = reshape(lambda_Init, lambdaDim, nStages);

% eta_Init
F_FuncObj_map = self.OCPEC.FuncObj.F.map(nStages);
ETA_Init = full(F_FuncObj_map(X_Init, U_Init, LAMBDA_Init));

%% group and reshape
z_Init = reshape([X_Init; U_Init; LAMBDA_Init; ETA_Init], [], 1);

end