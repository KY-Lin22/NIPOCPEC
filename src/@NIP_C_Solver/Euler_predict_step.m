function [gamma_h_Eu_j, gamma_c_Eu_j, gamma_g_Eu_j, z_Eu_j, Info] = Euler_predict_step(self,...
    gamma_h, gamma_c, gamma_g, z, p, p_j)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%
NLP = self.NLP;
Y_Node = cumsum([NLP.Dim.h, NLP.Dim.c, NLP.Dim.g, NLP.Dim.z]);

%% Function, Jacobian, Hessian, KKT matrix, and Sensitivity matrix evaluation
timeStart_funcGradEval = tic;

% Function and gradient evaluation
c = full(NLP.FuncObj.c(z));
g = full(NLP.FuncObj.g(z, p(1)));
h_grad = sparse(NLP.FuncObj.h_grad(z));
c_grad = sparse(NLP.FuncObj.c_grad(z));
g_grad = sparse(NLP.FuncObj.g_grad(z, p(1)));
% hessian
LAG_hessian = sparse(NLP.FuncObj.LAG_hessian(z, p(1)));
% smooth FB gradient
PSI_c_grad_dual  = sparse(self.FuncObj.PSI_c_grad_dual(gamma_c, c, p(2)));
PSI_c_grad_ineq  = sparse(self.FuncObj.PSI_c_grad_ineq(gamma_c, c, p(2)));
PSI_c_grad_sigma = sparse(self.FuncObj.PSI_c_grad_sigma(gamma_c, c, p(2)));
PSI_g_grad_dual  = sparse(self.FuncObj.PSI_g_grad_dual(gamma_g, g, p(2)));
PSI_g_grad_ineq  = sparse(self.FuncObj.PSI_g_grad_ineq(gamma_g, g, p(2)));
PSI_g_grad_sigma = sparse(self.FuncObj.PSI_g_grad_sigma(gamma_g, g, p(2)));
% KKT matrix
KKT_Matrix = self.evaluateKKT_Matrix(h_grad, c_grad, g_grad, LAG_hessian,...
    PSI_c_grad_dual, PSI_c_grad_ineq, PSI_g_grad_dual, PSI_g_grad_ineq);
% sensitivity matrix
Sensitivity_Matrix = self.evaluateSensitivity_Matrix(PSI_c_grad_sigma,...
    PSI_g_grad_ineq, PSI_g_grad_sigma);

timeElasped_funcGradEval = toc(timeStart_funcGradEval);

%% step evaluation
timeStart_stepEval = tic;

% compute Y_Eu_j
Y = [gamma_h; gamma_c; gamma_g; z];
Y_Eu_j = Y - KKT_Matrix\(Sensitivity_Matrix * (p_j - p));

timeElasped_stepEval = toc(timeStart_stepEval);

%% create output
% extract primal and dual variable from Y_Eu_j
gamma_h_Eu_j = Y_Eu_j(            1 : Y_Node(1), :); 
gamma_c_Eu_j = Y_Eu_j(Y_Node(1) + 1 : Y_Node(2), :) ;
gamma_g_Eu_j = Y_Eu_j(Y_Node(2) + 1 : Y_Node(3), :);
z_Eu_j       = Y_Eu_j(Y_Node(3) + 1 : Y_Node(4), :);
% info
Info.timeElasped_funcGradEval = timeElasped_funcGradEval;
Info.timeElasped_stepEval = timeElasped_stepEval;

end

