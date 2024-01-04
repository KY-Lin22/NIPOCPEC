function [gamma_h_j, gamma_c_j, gamma_g_j, z_j, Info] = ...
    Newton_correct_step(self, gamma_h_Eu_j, gamma_c_Eu_j, gamma_g_Eu_j, z_Eu_j, p_j)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
%
NLP = self.NLP;
Option = self.Option;
Y_Node = cumsum([NLP.Dim.h, NLP.Dim.c, NLP.Dim.g, NLP.Dim.z]);
timeElasped_funcGradEval = 0;
timeElasped_stepEval = 0;

%% Newton Corrector Step
% initialization
i = 1; % Newton corrector step counter
gamma_h_Ne = gamma_h_Eu_j;
gamma_c_Ne = gamma_c_Eu_j;
gamma_g_Ne = gamma_g_Eu_j;
z_Ne       = z_Eu_j;

while true
    %% KKT residual and matrix evaluation
    timeStart_funcGradEval_i = tic;
    
    [KKT_Residual_Ne, KKT_Matrix_Ne] = evaluateKKT_residual_matrix_correct_step(self,...
        gamma_h_Ne, gamma_c_Ne, gamma_g_Ne, z_Ne, p_j);
    
    timeElasped_funcGradEval_i = toc(timeStart_funcGradEval_i);
    timeElasped_funcGradEval = timeElasped_funcGradEval + timeElasped_funcGradEval_i;
    
    %% step evaluation
    timeStart_stepEval_i = tic;
    
    dY_Ne = - KKT_Matrix_Ne\(KKT_Residual_Ne);
    dgamma_h_Ne = dY_Ne(            1 : Y_Node(1), :);
    dgamma_c_Ne = dY_Ne(Y_Node(1) + 1 : Y_Node(2), :);
    dgamma_g_Ne = dY_Ne(Y_Node(2) + 1 : Y_Node(3), :);
    dz_Ne       = dY_Ne(Y_Node(3) + 1 : Y_Node(4), :);
    
    % evaluate new iterate
    switch Option.NewtonCorrection.employLineSearch
        case 0
            % full step
            gamma_h_Ne_j = gamma_h_Ne + dgamma_h_Ne;
            gamma_c_Ne_j = gamma_c_Ne + dgamma_c_Ne;
            gamma_g_Ne_j = gamma_g_Ne + dgamma_g_Ne;
            z_Ne_j       = z_Ne       + dz_Ne;
            terminalStatus = 1;
        case 1
            % line search
            [gamma_h_Ne_j, gamma_c_Ne_j, gamma_g_Ne_j, z_Ne_j, lineSearchStatus] =...
                LineSearch_Merit_correct_step(self, KKT_Residual_Ne,...
                gamma_h_Ne, gamma_c_Ne, gamma_g_Ne, z_Ne,...
                dgamma_h_Ne, dgamma_c_Ne, dgamma_g_Ne, dz_Ne, p_j);
%             terminalStatus = lineSearchStatus;            
            if lineSearchStatus ~= 1
                % line search fails, use full step
                gamma_h_Ne_j = gamma_h_Ne + dgamma_h_Ne;
                gamma_c_Ne_j = gamma_c_Ne + dgamma_c_Ne;
                gamma_g_Ne_j = gamma_g_Ne + dgamma_g_Ne;
                z_Ne_j       = z_Ne       + dz_Ne;
            end
            terminalStatus = 1;
    end  
    timeElasped_stepEval_i = toc(timeStart_stepEval_i);
    timeElasped_stepEval = timeElasped_stepEval + timeElasped_stepEval_i;   
    
    %% check termination
    if i == Option.NewtonCorrection.StepNum
        break        
    else
        % prepare for next Newton correction step
        i = i + 1;   
        gamma_h_Ne = gamma_h_Ne_j;
        gamma_c_Ne = gamma_c_Ne_j;
        gamma_g_Ne = gamma_g_Ne_j;
        z_Ne       = z_Ne_j;
    end
end

%% create output
gamma_h_j = gamma_h_Ne_j; 
gamma_c_j = gamma_c_Ne_j;
gamma_g_j = gamma_g_Ne_j;
z_j       = z_Ne_j;
% info
Info.timeElasped_funcGradEval = timeElasped_funcGradEval;
Info.timeElasped_stepEval = timeElasped_stepEval;
Info.terminalStatus = terminalStatus;
end

%% sub function: evaluate KKT residual and matrix
function [KKT_Residual, KKT_Matrix] = evaluateKKT_residual_matrix_correct_step(self, gamma_h, gamma_c, gamma_g, z, p)

% Function and gradient evaluation
h = full(self.NLP.FuncObj.h(z));
c = full(self.NLP.FuncObj.c(z));
g = full(self.NLP.FuncObj.g(z, p(1)));
J_grad = full(self.NLP.FuncObj.J_grad(z));
h_grad = sparse(self.NLP.FuncObj.h_grad(z));
c_grad = sparse(self.NLP.FuncObj.c_grad(z));
g_grad = sparse(self.NLP.FuncObj.g_grad(z, p(1)));
% hessian
LAG_hessian = sparse(self.NLP.FuncObj.LAG_hessian(z, p(1)));
% smooth FB function and its gradient
PSI_c = full(self.FuncObj.PSI_c(gamma_c, c, p(2)));
PSI_g = full(self.FuncObj.PSI_g(gamma_g, g, p(2)));
PSI_c_grad_dual  = sparse(self.FuncObj.PSI_c_grad_dual(gamma_c, c, p(2)));
PSI_c_grad_ineq  = sparse(self.FuncObj.PSI_c_grad_ineq(gamma_c, c, p(2)));
PSI_g_grad_dual  = sparse(self.FuncObj.PSI_g_grad_dual(gamma_g, g, p(2)));
PSI_g_grad_ineq  = sparse(self.FuncObj.PSI_g_grad_ineq(gamma_g, g, p(2)));
% KKT residual and KKT matrix
LAG_grad_z = J_grad + gamma_h' * h_grad - gamma_c' * c_grad  - gamma_g' * g_grad;
KKT_Residual = [h; PSI_c; PSI_g; LAG_grad_z'];
KKT_Matrix = self.evaluateKKT_Matrix(h_grad, c_grad, g_grad, LAG_hessian,...
    PSI_c_grad_dual, PSI_c_grad_ineq, PSI_g_grad_dual, PSI_g_grad_ineq);
end

%% sub function: merit line search
function [gamma_h_j, gamma_c_j, gamma_g_j, z_j, lineSearchStatus] =...
    LineSearch_Merit_correct_step(self, KKT_Residual,...
    gamma_h, gamma_c, gamma_g, z, ...
    dgamma_h, dgamma_c, dgamma_g, dz, p)

% load parameter
stepSize_min = self.Option.NewtonCorrection.stepSize_Min;
stepSize_decayRate = self.Option.NewtonCorrection.stepSize_DecayRate;
nu_D = self.Option.NewtonCorrection.nu_D;
% merit at current iterate Y
merit = 1/2 * norm(KKT_Residual, 2);
% backtracking line search
has_found_new_iterate = false;
stepSize_init = 1;
while ~has_found_new_iterate
    %% Step 1: estimate trial stepsize, iterate, and merit
    stepSize_trial = max([stepSize_init, stepSize_min]);
    gamma_h_trial = gamma_h + stepSize_trial * dgamma_h;
    gamma_c_trial = gamma_c + stepSize_trial * dgamma_c;
    gamma_g_trial = gamma_g + stepSize_trial * dgamma_g;
    z_trial       = z       + stepSize_trial * dz;
    % function and Jacobian
    h_trial = full(self.NLP.FuncObj.h(z_trial));
    c_trial = full(self.NLP.FuncObj.c(z_trial));
    g_trial = full(self.NLP.FuncObj.g(z_trial, p(1)));
    J_grad_trial = full(self.NLP.FuncObj.J_grad(z_trial));
    h_grad_trial = sparse(self.NLP.FuncObj.h_grad(z_trial));
    c_grad_trial = sparse(self.NLP.FuncObj.c_grad(z_trial));
    g_grad_trial = sparse(self.NLP.FuncObj.g_grad(z_trial, p(1)));
    PSI_c_trial = full(self.FuncObj.PSI_c(gamma_c_trial, c_trial, p(2)));
    PSI_g_trial = full(self.FuncObj.PSI_g(gamma_g_trial, g_trial, p(2)));
    LAG_grad_z_trial = J_grad_trial + gamma_h_trial' * h_grad_trial...
        - gamma_c_trial' * c_grad_trial  - gamma_g_trial' * g_grad_trial;
    % KKT residual
    KKT_Residual_trial = [h_trial; PSI_c_trial; PSI_g_trial; LAG_grad_z_trial'];
    % merit
    merit_trial = 1/2 * norm(KKT_Residual_trial, 2);
    
    %% Step 2: check sufficient decrease condition
    if merit_trial <= (1 - 2 * stepSize_trial * nu_D) * merit
        has_found_new_iterate = true;
        status = 1;
    end
    %% Step 3: checking min stepsize
    if ~has_found_new_iterate
        if stepSize_trial == stepSize_min
            % linesearch fails on the min stepsize, break backtracking linesearch procedure
            status = 0;
            break
        else
            % estimate a smaller stepsize
            stepSize_init = stepSize_decayRate * stepSize_init;
        end
    end
end
% organize output
lineSearchStatus = status;
switch status
    case 0
        % fail, return the previous one
        gamma_h_j = gamma_h;
        gamma_c_j = gamma_c;
        gamma_g_j = gamma_g;
        z_j = z;
    case 1
        % success, return the new iterate
        gamma_h_j = gamma_h_trial;
        gamma_c_j = gamma_c_trial;
        gamma_g_j = gamma_g_trial;
        z_j = z_trial;        
end
end

