function FuncObj = createFuncObj(self, nlp)
%createFuncObj
%   create function object
% output FuncObj is a structure with fields:
%        J, h, c, g: cost and constraint function
%        J_grad, h_grad, c_grad, g_grad: cost and constraint Jacobian 
%        LAG_hessian: Lagrangian Hessian
import casadi.*
FuncObj = struct();
%% function
% cost function
FuncObj.J = Function('J', {nlp.z}, {nlp.J}, {'z'}, {'J'});

% constraint
FuncObj.h = Function('h', {nlp.z}, {nlp.h}, {'z'}, {'h'});
FuncObj.c = Function('c', {nlp.z}, {nlp.c}, {'z'}, {'c'});
FuncObj.g = Function('g', {nlp.z, nlp.s}, {nlp.g}, {'z', 's'}, {'g'});

%% Jacobian
% cost function
FuncObj.J_grad = Function('J_grad', {nlp.z}, {nlp.J_grad}, {'z'}, {'J_grad'});

% constraint
FuncObj.h_grad = Function('h_grad', {nlp.z}, {nlp.h_grad}, {'z'}, {'h_grad'});
FuncObj.c_grad = Function('c_grad', {nlp.z}, {nlp.c_grad}, {'z'}, {'c_grad'});
FuncObj.g_grad = Function('g_grad', {nlp.z, nlp.s}, {nlp.g_grad}, {'z', 's'}, {'g_grad'});

%% Hessian
% Gauss-Newton
FuncObj.LAG_hessian = Function('LAG_hessian', {nlp.z, nlp.s}, {nlp.LAG_hessian}, {'z', 's'}, {'LAG_hessian'});

end

