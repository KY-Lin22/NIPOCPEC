function [gamma_h_j, gamma_c_j, gamma_g_j, z_j, Info] = ...
    Newton_correct_step(self, gamma_h_Eu_j, gamma_c_Eu_j, gamma_g_Eu_j, z_Eu_j, p_j)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
%
Y_Node = cumsum([self.NLP.Dim.h, self.NLP.Dim.c, self.NLP.Dim.g, self.NLP.Dim.z]);
timeElasped_funcEval = 0;
timeElasped_stepEval = 0;

%% Newton Corrector Step
% initialization
gamma_h = gamma_h_Eu_j;
gamma_c = gamma_c_Eu_j;
gamma_g = gamma_g_Eu_j;
z       = z_Eu_j;

for i = 1 : self.Option.NewtonCorrection.StepNum
    %% KKT residual and matrix evaluation at current iterate
    timeStart_funcEval_i = tic;
    
    % Function and gradient evaluation
    h = full(self.NLP.FuncObj.h(z));
    c = full(self.NLP.FuncObj.c(z));
    g = full(self.NLP.FuncObj.g(z, p_j(1)));
    J_grad = full(self.NLP.FuncObj.J_grad(z));
    h_grad = sparse(self.NLP.FuncObj.h_grad(z));
    c_grad = sparse(self.NLP.FuncObj.c_grad(z));
    g_grad = sparse(self.NLP.FuncObj.g_grad(z, p_j(1)));
    % hessian
    LAG_hessian = sparse(self.NLP.FuncObj.LAG_hessian(z, p_j(1)));
    % smooth FB function and its gradient
    PSI_c = full(self.FuncObj.PSI_c(gamma_c, c, p_j(2)));
    PSI_g = full(self.FuncObj.PSI_g(gamma_g, g, p_j(2)));
    PSI_c_grad_dual  = sparse(self.FuncObj.PSI_c_grad_dual(gamma_c, c, p_j(2)));
    PSI_c_grad_ineq  = sparse(self.FuncObj.PSI_c_grad_ineq(gamma_c, c, p_j(2)));
    PSI_g_grad_dual  = sparse(self.FuncObj.PSI_g_grad_dual(gamma_g, g, p_j(2)));
    PSI_g_grad_ineq  = sparse(self.FuncObj.PSI_g_grad_ineq(gamma_g, g, p_j(2)));
    % KKT residual and KKT matrix
    LAG_grad_z = J_grad + gamma_h' * h_grad - gamma_c' * c_grad  - gamma_g' * g_grad;
    KKT_Residual = [h; PSI_c; PSI_g; LAG_grad_z'];
    KKT_Matrix = self.evaluateKKT_Matrix(h_grad, c_grad, g_grad, LAG_hessian,...
        PSI_c_grad_dual, PSI_c_grad_ineq, PSI_g_grad_dual, PSI_g_grad_ineq);
    
    timeElasped_funcEval_i = toc(timeStart_funcEval_i);
    timeElasped_funcEval = timeElasped_funcEval + timeElasped_funcEval_i;

    %% step evaluation
    timeStart_stepEval_i = tic;
    
    dY = - KKT_Matrix\(KKT_Residual);
    dgamma_h = dY(            1 : Y_Node(1), :);
    dgamma_c = dY(Y_Node(1) + 1 : Y_Node(2), :);
    dgamma_g = dY(Y_Node(2) + 1 : Y_Node(3), :);
    dz       = dY(Y_Node(3) + 1 : Y_Node(4), :);
    
    % evaluate new iterate
    gamma_h = gamma_h + dgamma_h;
    gamma_c = gamma_c + dgamma_c;
    gamma_g = gamma_g + dgamma_g;
    z       = z       + dz;
    
    timeElasped_stepEval_i = toc(timeStart_stepEval_i);
    timeElasped_stepEval = timeElasped_stepEval + timeElasped_stepEval_i; 

end

%% create output
gamma_h_j = gamma_h; 
gamma_c_j = gamma_c;
gamma_g_j = gamma_g;
z_j       = z;

timeStart_funcEval_info = tic;
% evaluate function and gradient
J_j = full(self.NLP.FuncObj.J(z_j));
h_j = full(self.NLP.FuncObj.h(z_j));
c_j = full(self.NLP.FuncObj.c(z_j));
g_j = full(self.NLP.FuncObj.g(z_j, p_j(1)));
J_grad_j = full(self.NLP.FuncObj.J_grad(z_j));
h_grad_j = sparse(self.NLP.FuncObj.h_grad(z_j));
c_grad_j = sparse(self.NLP.FuncObj.c_grad(z_j));
g_grad_j = sparse(self.NLP.FuncObj.g_grad(z_j, p_j(1)));
% evaluate KKT error
LAG_grad_z_j = J_grad_j + gamma_h_j' * h_grad_j - gamma_c_j' * c_grad_j  - gamma_g_j' * g_grad_j;
[KKT_error_primal_j, KKT_error_dual_j, KKT_error_dual_scaled_j,...
    KKT_error_complementary_j, KKT_error_complementary_scaled_j, KKT_error_total_j] = ...
    self.evaluateKKT_error(gamma_h_j, gamma_c_j, gamma_g_j, h_j, c_j, g_j, LAG_grad_z_j);
% evaluate VI natural residual
VI_natural_residual_j  = self.evaluateNaturalResidual(z_j);

timeElasped_funcEval_info = toc(timeStart_funcEval_info);
timeElasped_funcEval = timeElasped_funcEval + timeElasped_funcEval_info;

% info: quantities
Info.cost = J_j;
Info.KKT_error_primal               = KKT_error_primal_j;
Info.KKT_error_dual                 = KKT_error_dual_j;
Info.KKT_error_dual_scaled          = KKT_error_dual_scaled_j;
Info.KKT_error_complementary        = KKT_error_complementary_j;
Info.KKT_error_complementary_scaled = KKT_error_complementary_scaled_j;
Info.KKT_error_total                = KKT_error_total_j;
Info.VI_natural_residual            = VI_natural_residual_j;
% info: time
Info.time.funcEval = timeElasped_funcEval;
Info.time.stepEval = timeElasped_stepEval;
end

