function FuncObj = createFuncObj(self)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*
FuncObj = struct();

%% initialize other MX-type symbolic variable and parameter
% dual variable
gamma_h = MX.sym('gamma_h', self.NLP.Dim.h, 1);
gamma_c = MX.sym('gamma_c', self.NLP.Dim.c, 1);
gamma_g = MX.sym('gamma_g', self.NLP.Dim.g, 1);
% smooth FB parameter
sigma = MX.sym('sigma', 1, 1);
% vector that collects variable and parameter
Y = [gamma_h; gamma_c; gamma_g; self.NLP.z];
p = [self.NLP.s; sigma];
% dY and dp
dgamma_h = MX.sym('dgamma_h', self.NLP.Dim.h, 1);
dgamma_c = MX.sym('dgamma_c', self.NLP.Dim.c, 1);
dgamma_g = MX.sym('dgamma_g', self.NLP.Dim.g, 1);
dz = MX.sym('dz', self.NLP.Dim.z, 1);
dY = [dgamma_h; dgamma_c; dgamma_g; dz];
dp = MX.sym('dp', 2, 1);

%% NLP function, Jacobian, and Hessian
% function
FuncObj.J = Function('J', {self.NLP.z}, {self.NLP.J}, {'z'}, {'J'});
FuncObj.h = Function('h', {self.NLP.z}, {self.NLP.h}, {'z'}, {'h'});
FuncObj.c = Function('c', {self.NLP.z}, {self.NLP.c}, {'z'}, {'c'});
FuncObj.g = Function('g', {self.NLP.z, self.NLP.s}, {self.NLP.g}, {'z', 's'}, {'g'});
% Jacobian
J_grad = jacobian(self.NLP.J, self.NLP.z);
h_grad = jacobian(self.NLP.h, self.NLP.z);
c_grad = jacobian(self.NLP.c, self.NLP.z);
g_grad = jacobian(self.NLP.g, self.NLP.z);
% Hessian
[J_hessian, ~] = hessian(self.NLP.J, self.NLP.z);

%% NLP Lagrangian
% function
LAG = self.NLP.J + gamma_h' * self.NLP.h - gamma_c' * self.NLP.c - gamma_g' * self.NLP.g;
% jacobian
LAG_grad = jacobian(LAG, self.NLP.z);
FuncObj.LAG_grad = Function('LAG_grad', {Y, p}, {LAG_grad}, {'Y', 'p'}, {'LAG_grad'});
% Hessian
switch self.Option.KKT.Hessian_approximation
    case 'Exact'
        % exact Hessian        
        [LAG_hessian, ~] = hessian(LAG, self.NLP.z);
    case 'Gauss_Newton'
        % Gauss-Newton Hessian approximation
        LAG_hessian = J_hessian;
    otherwise
        error('specified Hessian approximation method is not supported')
end

%% perturbed system of equation for complementarity between inequality constraint and its dual variable 
% create FuncObj of the smooth FB function 
a_SX = SX.sym('a', 1, 1);
b_SX = SX.sym('b', 1, 1);
sigma_SX = SX.sym('sigma', 1, 1);
psi = sqrt(a_SX^2 + b_SX^2 + sigma_SX^2) - a_SX - b_SX;
psi_FuncObj = Function('psi', {a_SX, b_SX, sigma_SX}, {psi}, {'a', 'b', 'sigma'}, {'psi'});

% map complementarity condition in KKT into perturbed system of equation
psi_FuncObj_map_c = psi_FuncObj.map(self.NLP.Dim.c);
PSI_c = (psi_FuncObj_map_c(gamma_c', self.NLP.c', sigma))';
PSI_c_grad_gamma_c = jacobian(PSI_c, gamma_c);
PSI_c_grad_z = jacobian(PSI_c, self.NLP.z);
PSI_c_grad_sigma = jacobian(PSI_c, sigma);

psi_FuncObj_map_g = psi_FuncObj.map(self.NLP.Dim.g);
PSI_g = (psi_FuncObj_map_g(gamma_g', self.NLP.g', sigma))';
PSI_g_grad_gamma_g = jacobian(PSI_g, gamma_g);
PSI_g_grad_z = jacobian(PSI_g, self.NLP.z);
PSI_g_grad_s = jacobian(PSI_g, self.NLP.s);
PSI_g_grad_sigma = jacobian(PSI_g, sigma);

%% line search in non-interior-point method
% constraint violation M
M = [self.NLP.h; PSI_c; PSI_g];
FuncObj.M = Function('M', {Y, p}, {M}, {'Y', 'p'}, {'M'});
% directional derivative of cost 
J_grad_times_dz = J_grad * dz;
FuncObj.J_grad_times_dz = Function('J_grad_times_dz', {Y, dY}, {J_grad_times_dz}, {'Y', 'dY'}, {'J_grad_times_dz'});

%% KKT residual and matrix
% KKT residual
KKT_residual = [self.NLP.h; PSI_c; PSI_g; LAG_grad'];
FuncObj.KKT_residual = Function('KKT_residual', {Y, p}, {KKT_residual}, {'Y', 'p'}, {'KKT_residual'});
% KKT matrix
nu_h = self.Option.RegParam.nu_h;
nu_c = self.Option.RegParam.nu_c;
nu_g = self.Option.RegParam.nu_g;
nu_H = self.Option.RegParam.nu_H;

KKT_matrix = ...
    [-nu_h * MX.eye(self.NLP.Dim.h),     MX(self.NLP.Dim.h, self.NLP.Dim.c),                  MX(self.NLP.Dim.h, self.NLP.Dim.g),                 h_grad;...
    MX(self.NLP.Dim.c, self.NLP.Dim.h),  PSI_c_grad_gamma_c - nu_c * MX.eye(self.NLP.Dim.c),  MX(self.NLP.Dim.c, self.NLP.Dim.g),                 PSI_c_grad_z;...
    MX(self.NLP.Dim.g, self.NLP.Dim.h),  MX(self.NLP.Dim.g, self.NLP.Dim.c),                  PSI_g_grad_gamma_g - nu_g * MX.eye(self.NLP.Dim.g), PSI_g_grad_z;...
    h_grad',                             -c_grad',                                            -g_grad',                                           LAG_hessian + nu_H * MX.eye(self.NLP.Dim.z)];

%% sensitivity matrix (w.r.t. parameter)
% formulate sensitivity matrix
sensitivity_matrix = ...
    [MX(self.NLP.Dim.h, 1), MX(self.NLP.Dim.h, 1);...
    MX(self.NLP.Dim.c, 1),  PSI_c_grad_sigma;...
    PSI_g_grad_s,           PSI_g_grad_sigma;...
    MX(self.NLP.Dim.z, 1),  MX(self.NLP.Dim.z, 1)];

%% Euler step and Newton step
dY_Euler = -KKT_matrix\(KKT_residual + sensitivity_matrix * dp);
dY_Newton = -KKT_matrix\KKT_residual;

FuncObj.dY_Euler = Function('dY_Euler', {Y, p, dp}, {dY_Euler}, {'Y', 'p', 'dp'}, {'dY_Euler'});
FuncObj.dY_Newton = Function('dY_Newton', {Y, p}, {dY_Newton}, {'Y', 'p'}, {'dY_Newton'});

end