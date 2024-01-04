function nlp = createNLP_box_constraint(self, OCPEC)
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
%        (3) J: cost function J = sum(J_n) and J_n = J_stage_n
%            with J_stage_n:   stage (and terminal) cost defined in OCPEC
%        (4) h: equality constraint arranged in a stagewise manner
%            h = [h_1;...h_n;...h_N] and h_n = [f_n; F_n - eta_n; C_n]   
%            with f_n:   discretized state equation f defined in OCPEC
%                 F_n:   VI function
%                 C_n:   path equality constraints C defined in OCPEC   
%        (5) c: inequality constraint (without s) arranged in a stagewise manner
%            c = [c_1;...c_n;...c_N] and c_n = [lambda_n - bl; bu - lambda_n; G_n] 
%            with G_n: path inequality constraints G defined in OCPEC
%                 bl, bu: VI set K, box constraint structure with lower bound bl and upper bound bu
%        (6) g: inequality constraint (with s) arranged in a stagewise manner
%            g = [g_1;...g_n;...g_N] and g_n = [s - (lambda_n - bl) * eta_n; s + (bu - lambda_n) * eta_n] 
%
% output: nlp is a structure with fields:
%         z: variable
%         s: parameter
%         J, h, c, g: cost and constraint function,                  
%         J_grad, h_grad, c_grad, g_grad: cost and constraint Jacobian
%         LAG_hessian: Lagrangian Hessian (Gauss-Newton)
%         Dim: problem size
%
import casadi.*

%% initialize NLP variable and function (stagewise, capital)
% define auxiliary variable dimension
eta_Dim = OCPEC.Dim.lambda;

% initialize problem variable 
X = SX.sym('X', OCPEC.Dim.x, OCPEC.nStages); 
U = SX.sym('U', OCPEC.Dim.u, OCPEC.nStages);
LAMBDA = SX.sym('LAMBDA', OCPEC.Dim.lambda, OCPEC.nStages);
ETA = SX.sym('ETA', eta_Dim, OCPEC.nStages);

Z = [X; U; LAMBDA; ETA]; % NLP variable in stagewise manner

% initialize problem parameter
s = SX.sym('s', 1, 1);

% initialize problem cost and constraint function
J = 0;
H = SX(OCPEC.Dim.x + eta_Dim + OCPEC.Dim.C, OCPEC.nStages); 
C = SX(2 * OCPEC.Dim.lambda + OCPEC.Dim.G, OCPEC.nStages);
G = SX(2 * OCPEC.Dim.lambda, OCPEC.nStages);

%% formulate NLP function
for n = 1 : OCPEC.nStages
    % load variable
    if n == 1
        x_nPrev = OCPEC.x0;
    else
        x_nPrev = X(:, n - 1);
    end    
    x_n = X(:, n);
    u_n = U(:, n);
    lambda_n = LAMBDA(:, n);      
    eta_n = ETA(:, n);      
    
    % NLP stage cost function
    J_stage_n = OCPEC.timeStep * OCPEC.FuncObj.L_S(x_n, u_n, lambda_n);
    if (n == OCPEC.nStages) && (size(OCPEC.L_T, 1) == 1)
        J_stage_n = J_stage_n + OCPEC.FuncObj.L_T(x_n);
    end    
    J = J + J_stage_n;    
    
    % discretized state equation f defined in OCPEC (implicit Euler method)
    f_n = x_nPrev - x_n + OCPEC.timeStep * OCPEC.FuncObj.f(x_n, u_n, lambda_n);  
    
    % summarize NLP constraint
    h_n = [f_n;...
        OCPEC.FuncObj.F(x_n, u_n, lambda_n) - eta_n;... % VI function
        OCPEC.FuncObj.C(x_n, u_n)...% path equality constraints C defined in OCPEC
        ];
    H(:, n) = h_n;      
    c_n = [lambda_n - OCPEC.bl;...
        OCPEC.bu - lambda_n;...
        OCPEC.FuncObj.G(x_n, u_n)...% path inequality constraints G defined in OCPEC
        ];
    C(:, n) = c_n;   
    g_n = [s * ones(OCPEC.Dim.lambda, 1) - (lambda_n - OCPEC.bl) .* eta_n;...
        s * ones(OCPEC.Dim.lambda, 1) + (OCPEC.bu - lambda_n) .* eta_n];
    G(:, n) = g_n;
end

%% formulate NLP variable and function (column, lowercase)
% problem variable
z = reshape(Z, (OCPEC.Dim.x + OCPEC.Dim.u + OCPEC.Dim.lambda + eta_Dim) * OCPEC.nStages, 1);

% constraint
h = reshape(H, (OCPEC.Dim.x + eta_Dim + OCPEC.Dim.C) * OCPEC.nStages, 1);
c = reshape(C, (2 * OCPEC.Dim.lambda + OCPEC.Dim.G) * OCPEC.nStages, 1);
g = reshape(G, (2 * OCPEC.Dim.lambda) * OCPEC.nStages, 1);

% size and node point of variable and function (NLP)
Dim.z = size(z, 1);
Dim.h = size(h, 1);
Dim.c = size(c, 1);
Dim.g = size(g, 1);
Dim.z_Node = cumsum([OCPEC.Dim.x, OCPEC.Dim.u, OCPEC.Dim.lambda, eta_Dim]); % node point after reshaping z into stagewise Z
Dim.h_Node = cumsum([OCPEC.Dim.x, eta_Dim, OCPEC.Dim.C]); % node point after reshaping h into stagewise H
Dim.c_Node = cumsum([OCPEC.Dim.lambda, OCPEC.Dim.lambda, OCPEC.Dim.G]); % node point after reshaping c into stagewise C
Dim.g_Node = cumsum([OCPEC.Dim.lambda, OCPEC.Dim.lambda]); % node point after reshaping g into stagewise G

%% formulate NLP Jacobian and Lagrangian Hessian
% cost function Jacobian
J_grad = jacobian(J, z); 
% constraint Jacobian
h_grad = jacobian(h, z);
c_grad = jacobian(c, z);
g_grad = jacobian(g, z);
% Gauss-Newton Hessian approximation
[J_hessian, ~] = hessian(J, z);
LAG_hessian = J_hessian;

%% create output struct
nlp = struct(...
    'z', z, 's', s,...
    'J', J, 'h', h, 'c', c, 'g', g,...
    'J_grad', J_grad, 'h_grad', h_grad, 'c_grad', c_grad, 'g_grad', g_grad,...
    'LAG_hessian', LAG_hessian,...
    'Dim', Dim);
end

