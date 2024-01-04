function FuncObj = createFuncObj(self)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*
NLP = self.NLP;

%% create FuncObj of the smooth FB function 
a = SX.sym('a', 1, 1);
b = SX.sym('b', 1, 1);
sigma = SX.sym('sigma', 1, 1);
psi = sqrt(a^2 + b^2 + sigma^2) - a - b;
psi_FuncObj = Function('psi', {a, b, sigma}, {psi}, {'a', 'b', 'sigma'}, {'psi'});

%% formulate smooth FB function for inequality constraint c, g
% initialize dual variable and inequality constraint
gamma_c = SX.sym('gamma_c', NLP.Dim.c, 1);
gamma_g = SX.sym('gamma_g', NLP.Dim.g, 1);
c = SX.sym('c', NLP.Dim.c, 1);
g = SX.sym('g', NLP.Dim.g, 1);

psi_FuncObj_map_c = psi_FuncObj.map(NLP.Dim.c);
PSI_c = (psi_FuncObj_map_c(gamma_c', c', sigma*ones(1, NLP.Dim.c)))';

psi_FuncObj_map_g = psi_FuncObj.map(NLP.Dim.g);
PSI_g = (psi_FuncObj_map_g(gamma_g', g', sigma*ones(1, NLP.Dim.g)))';

%% formulate smooth FB jacobian for inequality constraint c, g
PSI_c_grad_dual = jacobian(PSI_c, gamma_c);
PSI_c_grad_ineq = jacobian(PSI_c, c);
PSI_c_grad_sigma = jacobian(PSI_c, sigma);

PSI_g_grad_dual = jacobian(PSI_g, gamma_g);
PSI_g_grad_ineq = jacobian(PSI_g, g);
PSI_g_grad_sigma = jacobian(PSI_g, sigma);

%% create function object
FuncObj.PSI_c = Function('PSI_c', {gamma_c, c, sigma}, {PSI_c}, {'gamma_c', 'c', 'sigma'}, {'PSI_c'});
FuncObj.PSI_g = Function('PSI_g', {gamma_g, g, sigma}, {PSI_g}, {'gamma_g', 'g', 'sigma'}, {'PSI_g'});

FuncObj.PSI_c_grad_dual = Function('PSI_c_grad_dual', {gamma_c, c, sigma}, {PSI_c_grad_dual}, {'gamma_c', 'c', 'sigma'}, {'PSI_c_grad_dual'});
FuncObj.PSI_c_grad_ineq = Function('PSI_c_grad_ineq', {gamma_c, c, sigma}, {PSI_c_grad_ineq}, {'gamma_c', 'c', 'sigma'}, {'PSI_c_grad_ineq'});
FuncObj.PSI_c_grad_sigma = Function('PSI_c_grad_sigma', {gamma_c, c, sigma}, {PSI_c_grad_sigma}, {'gamma_c', 'c', 'sigma'}, {'PSI_c_grad_sigma'});

FuncObj.PSI_g_grad_dual = Function('PSI_g_grad_dual', {gamma_g, g, sigma}, {PSI_g_grad_dual}, {'gamma_g', 'g', 'sigma'}, {'PSI_g_grad_dual'});
FuncObj.PSI_g_grad_ineq = Function('PSI_g_grad_ineq', {gamma_g, g, sigma}, {PSI_g_grad_ineq}, {'gamma_g', 'g', 'sigma'}, {'PSI_g_grad_ineq'});
FuncObj.PSI_g_grad_sigma = Function('PSI_g_grad_sigma', {gamma_g, g, sigma}, {PSI_g_grad_sigma}, {'gamma_g', 'g', 'sigma'}, {'PSI_g_grad_sigma'});
end

