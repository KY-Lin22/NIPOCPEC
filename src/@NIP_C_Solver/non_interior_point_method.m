function [z_Opt, Info] = non_interior_point_method(self, z_Init, p)
%Stage 1 in NIP_C: solve NLP with given z_Init and p by non-interior-point method
% NLP has the form:
%  min  J(z),
%  s.t. h(z) = 0,
%       c(z) >= 0,
%       g(z, s) >= 0
% where: z is the variable,
%        s is the parameter,
%        J is the cost, and h, c, g are the constraint
% Syntax:
%          [z_Opt, Info] = non_interior_point_method(self, z_Init, p)
%          [z_Opt, Info] = self.non_interior_point_method(z_Init, p)
% Argument:
%          z_Init: double, NLP.Dim.z X 1, initial guess
%          p: double, given parameter, p = [s; sigma]
% Output:
%          z_Opt: double, NLP.Dim.z X 1, optimal solution found by solver
%          Info: struct, record the iteration information
%
%% iteration routine (z: previous iterate z_{k-1}, z_k: current iterate z_{k}) 
% time record
Time = struct('gradEval', 0, 'KKTEval', 0, 'searchDirection', 0, 'lineSearch', 0, 'else', 0, 'total', 0);
% Y node
Y_Node = cumsum([self.NLP.Dim.h, self.NLP.Dim.c, self.NLP.Dim.g, self.NLP.Dim.z]);
% extract parameter
s = p(1);
sigma = p(2);
% counter, beta, multipler, z, cost and constraint function
k = 1;
beta = self.Option.LineSearch.betaInit;
gamma_h = ones(self.NLP.Dim.h, 1);
gamma_c = ones(self.NLP.Dim.c, 1);
gamma_g = ones(self.NLP.Dim.g, 1);
z = z_Init;
J = full(self.NLP.FuncObj.J(z));
h = full(self.NLP.FuncObj.h(z));
c = full(self.NLP.FuncObj.c(z));
g = full(self.NLP.FuncObj.g(z, s));
% iteration loop
while true
    %% step 0: check iteration counter
    if k > self.Option.maxIterNum
        % failure case 1: exceed the maximum number of iteration
        terminalStatus = 0;
        terminalMsg = ['- Solver fails: ', 'because the maximum number of iteration exceeded'];
        break
    end    
    timeStart_total = tic;
    
    %% step 1: Jacobian and Hessian of previous iterate z
    timeStart_gradEval = tic;
    % cost and constraint Jacobian
    J_grad = full(self.NLP.FuncObj.J_grad(z));
    h_grad = sparse(self.NLP.FuncObj.h_grad(z));
    c_grad = sparse(self.NLP.FuncObj.c_grad(z));
    g_grad = sparse(self.NLP.FuncObj.g_grad(z, s));
    % Lagrangian Hessian
    LAG_hessian = sparse(self.NLP.FuncObj.LAG_hessian(z, s));
    % smooth FB function and its gradient
    PSI_c = full(self.FuncObj.PSI_c(gamma_c, c, sigma));
    PSI_g = full(self.FuncObj.PSI_g(gamma_g, g, sigma));
    PSI_c_grad_dual = sparse(self.FuncObj.PSI_c_grad_dual(gamma_c, c, sigma));
    PSI_c_grad_ineq = sparse(self.FuncObj.PSI_c_grad_ineq(gamma_c, c, sigma));
    PSI_g_grad_dual = sparse(self.FuncObj.PSI_g_grad_dual(gamma_g, g, sigma));
    PSI_g_grad_ineq = sparse(self.FuncObj.PSI_g_grad_ineq(gamma_g, g, sigma));   
    
    timeElasped_gradEval = toc(timeStart_gradEval);

    %% step 2: KKT residual, matrix, and error of previous iterate z
    timeStart_KKTEval = tic;
    % KKT residual
    LAG_grad_z = J_grad + gamma_h' * h_grad - gamma_c' * c_grad  - gamma_g' * g_grad;
    KKT_Residual = [h; PSI_c; PSI_g; LAG_grad_z'];   
    % KKT matrix
    KKT_Matrix = self.evaluateKKT_Matrix(h_grad, c_grad, g_grad, LAG_hessian,...
            PSI_c_grad_dual, PSI_c_grad_ineq, PSI_g_grad_dual, PSI_g_grad_ineq);
    % KKT error
    [KKT_error_primal, KKT_error_dual, KKT_error_dual_scaled,...
        KKT_error_complementary, KKT_error_complementary_scaled, KKT_error_total] = ...
        self.evaluateKKT_error(gamma_h, gamma_c, gamma_g, h, c, g, LAG_grad_z);

    timeElasped_KKTEval = toc(timeStart_KKTEval);
    
    %% step 3: search direction evaluation based on previous iterate z
    timeStart_SearchDirection = tic; 
    % solve sparse linear system using mldivide
    dY = KKT_Matrix\(-KKT_Residual); 
    % extract search direction
    dgamma_h = dY(            1 : Y_Node(1), 1);
    dgamma_c = dY(Y_Node(1) + 1 : Y_Node(2), 1);
    dgamma_g = dY(Y_Node(2) + 1 : Y_Node(3), 1);
    dz       = dY(Y_Node(3) + 1 : Y_Node(4), 1);
    % dYNorm (L_inf norm)
    dYNorm = norm(dY, inf); 
    timeElasped_searchDirection = toc(timeStart_SearchDirection);
    
    %% step 4: check whether we can terminate successfully based on the previous iterate z
    if KKT_error_total < self.Option.tol.KKT_error_total
        % Success case 1: the KKT error satisfies tolerance
        terminalStatus = 1;
        terminalMsg = ['- Solver succeeds: ', 'because the KKT error satisfies tolerance'];
        break
    elseif dYNorm < self.Option.tol.dYNorm
        % Success case 2: the norm of search direction satisfies tolerance
        terminalStatus = 1;
        terminalMsg = ['- Solver succeeds: ', 'because the norm of search direction satisfies tolerance'];   
        break
    elseif (KKT_error_primal <= self.Option.tol.KKT_error_primal) && ...
            (KKT_error_dual_scaled <= self.Option.tol.KKT_error_dual) && ...
            (KKT_error_complementary_scaled <= self.Option.tol.KKT_error_complementarity)
        % Success case 3: primal and dual error satisfy individual tolerance
        terminalStatus = 1;
        terminalMsg = ['- Solver succeeds: ', 'because primal, dual, and complementarity error satisfy individual tolerance']; 
        break        
    end
    
    %% step 5: merit line search
    [gamma_h_k, gamma_c_k, gamma_g_k, z_k, Info_LineSearch] = self.LineSearch_Merit(beta,...
            gamma_h, gamma_c, gamma_g, z, dgamma_h, dgamma_c, dgamma_g, dz,...
            s, sigma, J, h, PSI_c, PSI_g, J_grad);
    % check status
    if Info_LineSearch.status == 0
        % failure case 2: line search fails
        terminalStatus = -1;
        terminalMsg = ['- Solver fails: ', 'because merit line search reaches the min stepsize'];        
        break
    else
        % extract quantities (J, h, c, g) associated with z_k
        J_k = Info_LineSearch.J;
        h_k = Info_LineSearch.h;
        c_k = Info_LineSearch.c;
        g_k = Info_LineSearch.g;
        % line search quantities
        beta_k = Info_LineSearch.beta;
        stepSize = Info_LineSearch.stepSize;
        merit    = Info_LineSearch.merit;
    end    
    timeElasped_lineSearch = Info_LineSearch.time;
    
    %% step 6: record and print information of this iteration k
    timeElasped_total = toc(timeStart_total);
    Time.gradEval = Time.gradEval + timeElasped_gradEval;
    Time.KKTEval = Time.KKTEval + timeElasped_KKTEval;
    Time.searchDirection = Time.searchDirection + timeElasped_searchDirection;
    Time.lineSearch = Time.lineSearch + timeElasped_lineSearch;
    Time.total = Time.total + timeElasped_total;    
    % print
    if self.Option.printLevel == 2
        % head
        if mod(k, 10) == 1
            disp('----------------------------------------------------------------------------------------------------------------------------------------------------')
            headMsg = ' Iter |   cost   |  KKT(P)  |  KKT(D)  | KKT(D,s) |  KKT(C)  | KKT(C,s) | PSI(max) |  dYNorm  |   beta   | stepsize |  merit   | merit(t) | time(ms) |';
            disp(headMsg)
        end
        % previous iterate message
        prevIterMsg = ['  ',...
            num2str(k,'%10.3d'),' | ',...
            num2str(J,'%10.2e'), ' | ',...
            num2str(KKT_error_primal, '%10.2e'), ' | ',...
            num2str(KKT_error_dual, '%10.2e'), ' | ',...
            num2str(KKT_error_dual_scaled, '%10.2e'), ' | ',...
            num2str(KKT_error_complementary, '%10.2e'), ' | ',...
            num2str(KKT_error_complementary_scaled, '%10.2e'), ' | ',...
            num2str(max([norm(PSI_c, inf), norm(PSI_g, inf)]), '%10.2e'), ' | ',...
            num2str(dYNorm,'%10.2e'), ' | ',...
            num2str(beta_k,'%10.2e'), ' | ',...
            num2str(stepSize,'%10.2e'), ' | ',...
            num2str(merit(1),'%10.2e'), ' | ',...
            num2str(merit(2),'%10.2e'), ' | ',...
            num2str(1000 * timeElasped_total,'%10.2e'), ' | '];
        disp(prevIterMsg)
    end
    
    %% step 7: prepare next iteration
    k = k + 1;
    beta = beta_k;
    gamma_h = gamma_h_k;
    gamma_c = gamma_c_k;
    gamma_g = gamma_g_k;
    z = z_k;
    J = J_k;
    h = h_k;
    c = c_k; 
    g = g_k; 
end

%% return optimal solution and create information
% return previous iterate as solution
z_Opt = z;
% create Info (basic: time, iterNum, terminalStatus)
Time.else = Time.total - Time.searchDirection - Time.lineSearch - Time.KKTEval - Time.gradEval;
Info.Time = Time;
Info.iterNum = k - 1;
Info.terminalStatus = terminalStatus;
Info.terminalMsg = terminalMsg;
% create Info (corresponds to the solution z_Opt: dual variable, cost, KKT, natural residual)
Info.gamma_h = gamma_h;
Info.gamma_c = gamma_c;
Info.gamma_g = gamma_g;
Info.cost = J;
Info.KKT_error_primal               = KKT_error_primal;
Info.KKT_error_dual                 = KKT_error_dual;
Info.KKT_error_dual_scaled          = KKT_error_dual_scaled;
Info.KKT_error_complementary        = KKT_error_complementary;
Info.KKT_error_complementary_scaled = KKT_error_complementary_scaled;
Info.KKT_error_total                = KKT_error_total;
Info.VI_natural_residual            = self.evaluateNaturalResidual(z);
% display termination and solution message, then break rountie
if (self.Option.printLevel == 1) || (self.Option.printLevel == 2)
    disp('*--------------------------------------------- Solution Information ----------------------------------------------*')
    disp('1. Terminal Status')
    disp(Info.terminalMsg)
    disp('2. Iteration Process Message')
    disp(['- Iterations: .................................. ', num2str(Info.iterNum)])
    disp(['- TimeElapsed: ................................. ', num2str(Info.Time.total,'%10.3f'), 's'])
    disp(['- AverageTime: ................................. ', num2str(1000 * Info.Time.total /Info.iterNum,'%10.2f'), ' ms/Iter'])
    disp('3. Solution Message')
    disp(['- Cost: ........................................ ', num2str(Info.cost,'%10.3e'), '; '])
    disp(['- KKT (primal): ................................ ', num2str(Info.KKT_error_primal,'%10.3e'), '; '])
    disp(['- KKT (dual, unscaled / scaled): ............... ', num2str(Info.KKT_error_dual,'%10.3e'), ' / ', num2str(Info.KKT_error_dual_scaled,'%10.3e')  '; '])
    disp(['- KKT (complementary, unscaled / scaled): ...... ', num2str(Info.KKT_error_complementary,'%10.3e'), ' / ', num2str(Info.KKT_error_complementary_scaled,'%10.3e'), '; '])  
    disp(['- KKT (total): ................................. ', num2str(Info.KKT_error_total,'%10.3e'), '; '])
    disp(['- VI natural residual: ......................... ', num2str(Info.VI_natural_residual,'%10.3e'), '; '])
end

end

