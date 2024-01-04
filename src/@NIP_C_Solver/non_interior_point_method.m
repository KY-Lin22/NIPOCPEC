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
%% check option and input
% check option (TODO: write a function)
if self.Option.printLevel == 2
    self.Option.recordLevel = 1;
end
OCPEC = self.OCPEC;
NLP = self.NLP;
Option = self.Option;
% check input z_Init
if ~all(size(z_Init) == [NLP.Dim.z, 1])
    error('z_Init has wrong dimension')
end
% check parameter
if (p(1) < 0)
    error('relax parameter s (i.e, p_1) should be nonnegative')
end
if (p(2) < 0)
    error('FB smooth parameter sigma (i.e, p_2) should be nonnegative')
end

%% create record for time, log
% time
Time = struct('gradEval', 0, 'KKTEval', 0, 'searchDirection', 0, 'lineSearch', 0, 'else', 0, 'total', 0);
% log 
if Option.recordLevel == 1
    Log.cost      = zeros(Option.maxIterNum, 1); 
    Log.KKT_error = zeros(Option.maxIterNum, 6); % [primal, dual, dual_scaled, complementary, complementary_scaled, total] 
    Log.FB_max    = zeros(Option.maxIterNum, 2); % [PSI_c, PSI_g]
    Log.dYNorm    = zeros(Option.maxIterNum, 1);
    Log.beta      = zeros(Option.maxIterNum, 1); 
    Log.stepSize  = zeros(Option.maxIterNum, 1);
    Log.merit     = zeros(Option.maxIterNum, 2); % [merit, merit_k]
end

%% prepare the first iteration (z: previous iterate z_{k-1}, z_k: current iterate z_{k}) 
if (Option.printLevel == 1) || (Option.printLevel == 2)
    disp('**************************************************************************************************************************');
    disp('Initializing...');
end
% Y node
Y_Node = cumsum([NLP.Dim.h, NLP.Dim.c, NLP.Dim.g, NLP.Dim.z]);
% preproccess initial guess and extract parameter
z_Init = self.preprocessInitialGuess(z_Init);
s = p(1);
sigma = p(2);
% counter, beta, multipler, z, cost and constraint function
k = 1;
beta = Option.LineSearch.betaInit;
gamma_h = ones(NLP.Dim.h, 1);
gamma_c = ones(NLP.Dim.c, 1);
gamma_g = ones(NLP.Dim.g, 1);
z = z_Init;
J = full(NLP.FuncObj.J(z));
h = full(NLP.FuncObj.h(z));
c = full(NLP.FuncObj.c(z));
g = full(NLP.FuncObj.g(z, s));

%% iteration routine 
% begin iteration routine
if (Option.printLevel == 1) || (Option.printLevel == 2)
    disp('Computing the Optimal Solution...');
end
% iteration loop
while true
    %% step 0: check iteration counter
    if k > Option.maxIterNum
        % failure case 1: exceed the maximum number of iteration
        terminalStatus = 0;
        terminalMsg = ['- Solver fails: ', 'because the maximum number of iteration exceeded'];
        break
    end    
    timeStart_total = tic;
    
    %% step 1: Jacobian and Hessian of previous iterate z
    timeStart_gradEval = tic;
    % cost and constraint Jacobian
    J_grad = full(NLP.FuncObj.J_grad(z));
    h_grad = sparse(NLP.FuncObj.h_grad(z));
    c_grad = sparse(NLP.FuncObj.c_grad(z));
    g_grad = sparse(NLP.FuncObj.g_grad(z, s));
    % Lagrangian Hessian
    LAG_hessian = sparse(NLP.FuncObj.LAG_hessian(z, s));
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
    scaling_dual = max([Option.KKT_scaling_max,...
        norm([gamma_h; gamma_c; gamma_g], 1)/(NLP.Dim.h + NLP.Dim.c + NLP.Dim.g)])/Option.KKT_scaling_max;
    scaling_complementary = max([Option.KKT_scaling_max,...
        norm([gamma_c; gamma_g], 1)/(NLP.Dim.c + NLP.Dim.g)])/Option.KKT_scaling_max;
    KKT_error_primal = norm([h;...
        min([zeros(NLP.Dim.c, 1), c], [], 2);...
        min([zeros(NLP.Dim.g, 1), g], [], 2)], inf);
    KKT_error_dual = norm([LAG_grad_z';...
        min([zeros(NLP.Dim.c, 1), gamma_c], [], 2);...
        min([zeros(NLP.Dim.g, 1), gamma_g], [], 2)], inf); 
    KKT_error_dual_scaled = KKT_error_dual/scaling_dual;
    KKT_error_complementary = norm([c .* gamma_c; g .* gamma_g], inf);
    KKT_error_complementary_scaled = KKT_error_complementary/scaling_complementary;
    KKT_error_total = max([KKT_error_primal, KKT_error_dual_scaled, KKT_error_complementary_scaled]);

    timeElasped_KKTEval = toc(timeStart_KKTEval);
    
    %% step 3: search direction evaluation based on previous iterate z
    timeStart_SearchDirection = tic; 
    % solve linear system using mldivide(J is sparse)
    dY = KKT_Matrix\(-KKT_Residual); 
    % extract search direction
    dgamma_h = dY(            1 : Y_Node(1), 1);
    dgamma_c = dY(Y_Node(1) + 1 : Y_Node(2), 1);
    dgamma_g = dY(Y_Node(2) + 1 : Y_Node(3), 1);
    dz       = dY(Y_Node(3) + 1 : Y_Node(4), 1);
    % dYNorm (L2 norm scaled by nStages)
    dYNorm = norm(dY, 2)/OCPEC.nStages; 
    timeElasped_searchDirection = toc(timeStart_SearchDirection);
    
    %% step 4: check whether we can terminate successfully based on the previous iterate z
    if KKT_error_total < Option.tol.KKT_error_total
        % Success case 1: the KKT error satisfies tolerance
        terminalStatus = 1;
        terminalMsg = ['- Solver succeeds: ', 'because the KKT error satisfies tolerance'];
        break
    elseif dYNorm < Option.tol.dYNorm
        % Success case 2: the norm of search direction satisfies tolerance
        terminalStatus = 1;
        terminalMsg = ['- Solver succeeds: ', 'because the norm of search direction satisfies tolerance'];   
        break
    elseif (KKT_error_primal <= Option.tol.KKT_error_primal) && (KKT_error_dual_scaled <= Option.tol.KKT_error_dual)
        % Success case 3: primal and dual error satisfy individual tolerance
        terminalStatus = 1;
        terminalMsg = ['- Solver succeeds: ', 'because primal and dual error satisfy individual tolerance']; 
        break        
    end
    
    %% step 5: merit line search
    [gamma_h_k, gamma_c_k, gamma_g_k, z_k, Info_LineSearch] = self.LineSearch_Merit(beta,...
            gamma_h, gamma_c, gamma_g, z, dgamma_h, dgamma_c, dgamma_g, dz, p, J, h, PSI_c, PSI_g, J_grad);
    % check status
    if Info_LineSearch.status == 0
        % failure case 2: line search fails
        terminalStatus = -1;
        terminalMsg = ['- Solver fails: ', 'because merit line search reaches the min stepsize'];        
        break
    end    
    timeElasped_lineSearch = Info_LineSearch.time;
    % extract quantities (J, h, c, g) associated with z_k
    J_k = Info_LineSearch.J;
    h_k = Info_LineSearch.h;
    c_k = Info_LineSearch.c;
    g_k = Info_LineSearch.g;
    beta_k = Info_LineSearch.beta;
    stepSize = Info_LineSearch.stepSize;    
    merit    = Info_LineSearch.merit;
    
    %% step 6: record and print information of this iteration k
    timeElasped_total = toc(timeStart_total);
    Time.gradEval = Time.gradEval + timeElasped_gradEval;
    Time.KKTEval = Time.KKTEval + timeElasped_KKTEval;
    Time.searchDirection = Time.searchDirection + timeElasped_searchDirection;
    Time.lineSearch = Time.lineSearch + timeElasped_lineSearch;
    Time.total = Time.total + timeElasped_total;
    
    if Option.recordLevel == 1
        % record
        Log.cost(k)         = J;
        Log.KKT_error(k, :) = [KKT_error_primal, KKT_error_dual, KKT_error_dual_scaled, ...
            KKT_error_complementary, KKT_error_complementary_scaled, KKT_error_total];
        Log.FB_max(k, :)    = [norm(PSI_c, inf), norm(PSI_g, inf)];
        Log.dYNorm(k)       = dYNorm;
        Log.beta(k)         = beta_k;
        Log.stepSize(k)     = stepSize;        
        Log.merit(k, :)     = merit; 
        % print
        if Option.printLevel == 2
            % head
            if mod(k, 10) == 1
                disp('----------------------------------------------------------------------------------------------------------------------------------------------------')
                headMsg = ' Iter |   cost   |  KKT(P)  |  KKT(D)  | KKT(D,s) |  KKT(C)  | KKT(C,s) | PSI(max) |  dYNorm  |   beta   | stepsize |  merit   | merit(t) | time(ms) |';
                disp(headMsg)
            end 
            % previous iterate message
            prevIterMsg = ['  ',...
                num2str(k,'%10.3d'),' | ',...
                num2str(Log.cost(k),'%10.2e'), ' | ',...
                num2str(Log.KKT_error(k, 1), '%10.2e'), ' | ',...
                num2str(Log.KKT_error(k, 2), '%10.2e'), ' | ',...
                num2str(Log.KKT_error(k, 3), '%10.2e'), ' | ',...
                num2str(Log.KKT_error(k, 4), '%10.2e'), ' | ',...
                num2str(Log.KKT_error(k, 5), '%10.2e'), ' | ',...                
                num2str(max(Log.FB_max(k, :)), '%10.2e'), ' | ',...
                num2str(Log.dYNorm(k),'%10.2e'), ' | ',...
                num2str(Log.beta(k),'%10.2e'), ' | ',...
                num2str(Log.stepSize(k),'%10.2e'), ' | ',...
                num2str(Log.merit(k, 1),'%10.2e'), ' | ',...
                num2str(Log.merit(k, 2),'%10.2e'), ' | ',...
                num2str(1000 * timeElasped_total,'%10.2e'), ' | '];
            disp(prevIterMsg)            
        end
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
Info.cost    = J;
Info.KKT_error_primal               = KKT_error_primal;
Info.KKT_error_dual                 = KKT_error_dual;
Info.KKT_error_dual_scaled          = KKT_error_dual_scaled;
Info.KKT_error_complementary        = KKT_error_complementary;
Info.KKT_error_complementary_scaled = KKT_error_complementary_scaled;
Info.KKT_error_total                = KKT_error_total;
Info.VI_natural_residual            = self.evaluateNaturalResidual(z);
% create Info (log)
if Option.recordLevel == 1
    Info.Log.cost      = Log.cost(1 : k - 1, 1);
    Info.Log.KKT_error = Log.KKT_error(1 : k - 1, :);
    Info.Log.FB_max    = Log.FB_max(1 : k - 1, :);
    Info.Log.dYNorm    = Log.dYNorm(1 : k - 1, 1);
    Info.Log.beta      = Log.beta(1 : k - 1, 1);
    Info.Log.stepSize  = Log.stepSize(1 : k - 1, 1);
    Info.Log.merit     = Log.merit(1 : k - 1, :);
end
% display termination and solution message, then break rountie
if (Option.printLevel == 1) || (Option.printLevel == 2)
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
% end of the iteration routine 
if (Option.printLevel == 1) || (Option.printLevel == 2)
    disp('Done!')
    disp('**************************************************************************************************************************');
end

end

