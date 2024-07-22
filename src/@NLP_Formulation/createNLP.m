function nlp = createNLP(self, OCPEC)
%createNLP: formulate a NLP based on the given OCPEC (box constraint case)
%
% OCPEC has the form:
%  min  L_T(x) + int_0^T L_S(x, u, lambda) dt,
%  s.t. Dot{x} = f(x, u, lambda)
%       lambda \in SOL(K, F(x, u, lambda))
%       K := {lambda | bl <= lambda <= bu}
%       G(x, u) >= 0,
%       C(x, u) = 0,
%
% NLP has the form:
%  min  J(z),
%  s.t. h(z) = 0,
%       c(z) >= 0
%       g(z, s) >= 0
% where: (1) z: collects all the variables to be optimized and arranged in a stagewise manner
%            z = [z_1;...z_n;...z_N] and z_n = [x_n; u_n; lambda_n; eta_n] 
%            with x_n:      system state
%                 u_n:      control                 
%                 lambda_n: algebraic variable   
%                 eta_n:    auxiliary variable for VI function F       
%        (2) s: nonnegative relax parameter
%        (3) J: cost function J = sum(J_stage_n)*dt + J_terminal
%            with J_stage_n:   stage cost defined in OCPEC
%                 J_terminal:  terminal cost defined in OCPEC
%        (4) h: equality constraint arranged in a stagewise manner
%            h = [h_1;...h_n;...h_N] and h_n = [f_n; F_n - eta_n; C_n]   
%            with f_n:   discretized state equation f defined in OCPEC
%                 F_n:   VI function
%                 C_n:   path equality constraints C defined in OCPEC
%        (5) c: inequality constraint (without s) arranged in a stagewise manner
%            c = [c_1;...c_n;...c_N] and c_n = [equilibrium_constraint_without_param_n; G_n] 
%            with G_n: path inequality constraints G defined in OCPEC
%        (6) g: inequality constraint (with s) arranged in a stagewise manner
%            g = [g_1;...g_n;...g_N] and g_n = equilibrium_constraint_with_param_n 
% output: nlp is a structure with fields:
%         z: variable
%         s: parameter
%         J, h, c, g: cost and constraint function,                  
%         Dim: problem size
%
import casadi.*

%% initialize NLP variable (stagewise, capital)
% initialize problem variable (cell)
X_cell = cell(1, OCPEC.nStages);
U_cell = cell(1, OCPEC.nStages);
LAMBDA_cell = cell(1, OCPEC.nStages);
ETA_cell = cell(1, OCPEC.nStages);

for n = 1 : OCPEC.nStages
    X_cell{1, n} = MX.sym(['x_' num2str(n)], OCPEC.Dim.x, 1);
    U_cell{1, n} = MX.sym(['u_' num2str(n)], OCPEC.Dim.u, 1);
    LAMBDA_cell{1, n} = MX.sym(['lambda_' num2str(n)], OCPEC.Dim.lambda, 1);
    ETA_cell{1, n} = MX.sym(['eta_' num2str(n)], OCPEC.Dim.lambda, 1);
end

XPrev_cell = [OCPEC.x0, X_cell(1, 1 : end - 1)];
Z_cell = [X_cell; U_cell; LAMBDA_cell; ETA_cell];

% initialize problem variable (MX)
X = horzcat(X_cell{:});
U = horzcat(U_cell{:});
LAMBDA = horzcat(LAMBDA_cell{:});
ETA = horzcat(ETA_cell{:});
XPrev = horzcat(XPrev_cell{:});
Z = vertcat(Z_cell{:});

% initialize relaxation parameter
s = MX.sym('s', 1, 1);

%% mapping function object
% stage cost
L_S_map = OCPEC.FuncObj.L_S.map(OCPEC.nStages);
% DVI
f_map = OCPEC.FuncObj.f.map(OCPEC.nStages);
F_map = OCPEC.FuncObj.F.map(OCPEC.nStages);
% path inequality constraints G and equality constraint C defined in OCPEC
G_map = OCPEC.FuncObj.G.map(OCPEC.nStages);
C_map = OCPEC.FuncObj.C.map(OCPEC.nStages);
% relax equilibrium constraint
[constraint_without_param, constraint_with_param] = self.createRelaxedEquilibriumConstraint(OCPEC);
constraint_without_param_map = constraint_without_param.map(OCPEC.nStages);
constraint_with_param_map = constraint_with_param.map(OCPEC.nStages);

%% formulate NLP function (stagewise)
% cost
L_S_stage = L_S_map(X, U, LAMBDA);
L_T = OCPEC.FuncObj.L_T(X(:, end));
% DVI
f_stage = f_map(X, U, LAMBDA);
F_stage = F_map(X, U, LAMBDA);
% path inequality constraint G and equality constraint C
G_stage = G_map(X, U);
C_stage = C_map(X, U);
% relax equilibrium constraint
constraint_without_param_stage = constraint_without_param_map(LAMBDA, ETA);
constraint_with_param_stage = constraint_with_param_map(LAMBDA, ETA, s);

%% reshape NLP variable and function (column, lowercase)
% variable
z = reshape(Z, (OCPEC.Dim.x + OCPEC.Dim.u + 2 * OCPEC.Dim.lambda) * OCPEC.nStages, 1);
Dim.z_Node = cumsum([OCPEC.Dim.x, OCPEC.Dim.u, OCPEC.Dim.lambda, OCPEC.Dim.lambda]);
Dim.z = size(z, 1);

% cost function
J = sum(L_S_stage) * OCPEC.timeStep + L_T;

% equality constraint h = 0
h_stage = [...
    XPrev - X + f_stage * OCPEC.timeStep;...
    F_stage - ETA;...
    C_stage];
h = reshape(h_stage, (OCPEC.Dim.x + OCPEC.Dim.lambda + OCPEC.Dim.C) * OCPEC.nStages, 1);
Dim.h_Node = cumsum([OCPEC.Dim.x, OCPEC.Dim.lambda, OCPEC.Dim.C]);
Dim.h = size(h, 1);

% inequality constraint c >= 0
c_stage = [...
    constraint_without_param_stage;...
    G_stage];
c = reshape(c_stage, (size(constraint_without_param_stage, 1) + OCPEC.Dim.G) * OCPEC.nStages, 1);
Dim.c_Node = cumsum([size(constraint_without_param_stage, 1), OCPEC.Dim.G]);
Dim.c = size(c, 1);

% inequality constraint g >= 0
g_stage = constraint_with_param_stage;
g = reshape(g_stage, size(constraint_with_param_stage, 1) * OCPEC.nStages, 1);
Dim.g = size(g, 1);

%% create output struct
nlp = struct('z', z, 's', s,...
    'J', J, 'h', h, 'c', c, 'g', g,...
    'Dim', Dim);
end