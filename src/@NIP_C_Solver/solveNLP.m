function [z_Opt, Info] = solveNLP(self, z_Init, p_Init, p_End)
% solve NLP with given z_Init, p_Init, and p_End by non-interior-point continuation method
% NLP has the form:
%  min  J(z),
%  s.t. h(z) = 0,
%       c(z) >= 0,
%       g(z, s) >= 0
% where: z is the variable,
%        s is the parameter,
%        J is the cost, and h, c, g are the constraint
% Syntax:
%          [z_Opt, Info] = solveNLP(self, z_Init, p_Init, p_End)
%          [z_Opt, Info] = self.solveNLP(z_Init, p_Init, p_End)
% Argument:
%          z_Init: double, NLP.Dim.z X 1, initial guess
%          p_Init: double, problem parameter (initial) p = [s; sigma]
%          p_End: double, problem parameter (end) p = [s; sigma]
% Output:
%          z_Opt: double, NLP.Dim.z X 1, optimal solution found by solver
%          Info: struct, record the iteration information
%
%% check option and input
% check option (TODO)
% self.Option.printLevel = 0; % do not print stage 1 iteration information
%
NLP = self.NLP;
Option = self.Option;
% check input z_Init
if ~all(size(z_Init) == [NLP.Dim.z, 1])
    error('z_Init has wrong dimension')
end
% check parameter
if (p_Init(1) < 0) || (p_End(1) < 0)
    error('relax parameter s (i.e., p_1) should be nonnegative')
end
if (p_Init(2) < 0) || (p_End(2) < 0)
    error('FB smooth parameter sigma (i.e, p_2) should be nonnegative')
end
if p_Init(1) < p_End(1)
    error('s_Init should not smaller than s_End')
end
if p_Init(2) < p_End(2)
    error('sigma_Init should not smaller than sigma_End')
end
% load parameter
kappa_times = Option.Homotopy.kappa_times;
kappa_exp = Option.Homotopy.kappa_exp;

%% create record for time and log 
% evaluate the number of continuation step based on p_Init and p_End
p_test = p_Init;
continuationStepNum = 0;
while true
    if all(p_test == p_End)
        break
    else        
        p_trial = min([kappa_times .* p_test, p_test.^kappa_exp], [], 2);
        p_test = max([p_trial, p_End], [], 2);
        continuationStepNum = continuationStepNum + 1;
    end
end
% time 
Time = struct('non_interior_point', 0,...
    'Euler_predict_funcGradEval', 0, 'Euler_predict_stepEval', 0,...
    'Newton_correct_funcGradEval', 0, 'Newton_correct_stepEval',0,...
    'else', 0, 'total', 0);
% log
Log.param               = zeros(continuationStepNum + 1, 2); % [s, sigma]
Log.cost                = zeros(continuationStepNum + 1, 1);
Log.KKT_error           = zeros(continuationStepNum + 1, 6); % [primal, dual, dual_scaled, complementary, complementary_scaled, total]
Log.VI_natural_residual = zeros(continuationStepNum + 1, 1);
Log.timeElapsed         = zeros(continuationStepNum + 1, 1); % elapsed time in each continuation step

%% non-interior-point continuation method
while true
    %% stage 1: solve the first parameterized NLP by non-interior-point method
    % solve
    p = p_Init;
    [z, Info_Stage_1] = self.non_interior_point_method(z_Init, p);
    % extract quantities (dual variable, cost, KKT error, VI residual)
    gamma_h = Info_Stage_1.gamma_h;
    gamma_c = Info_Stage_1.gamma_c;
    gamma_g = Info_Stage_1.gamma_g;   
    J = Info_Stage_1.cost;
    KKT_error_primal               = Info_Stage_1.KKT_error_primal;
    KKT_error_dual                 = Info_Stage_1.KKT_error_dual;
    KKT_error_dual_scaled          = Info_Stage_1.KKT_error_dual_scaled;
    KKT_error_complementary        = Info_Stage_1.KKT_error_complementary;
    KKT_error_complementary_scaled = Info_Stage_1.KKT_error_complementary_scaled;
    KKT_error_total                = Info_Stage_1.KKT_error_total;
    VI_natural_residual            = Info_Stage_1.VI_natural_residual;
    % record  
    Time.non_interior_point = Info_Stage_1.Time.total;
    Time.total = Time.total + Info_Stage_1.Time.total;
    Log.param(1, :) = p';
    Log.cost(1, :)  = J;
    Log.KKT_error(1, :) = [KKT_error_primal, KKT_error_dual, KKT_error_dual_scaled,...
        KKT_error_complementary, KKT_error_complementary_scaled, KKT_error_total];
    Log.VI_natural_residual(1, :) = VI_natural_residual;
    Log.timeElapsed(1, :) = Info_Stage_1.Time.total;
    % print
    disp('---------------------------------------------------------------------------------------------------------------------------------------------')
    headMsg = ' StepNum |    s     |   sigma  |   cost   |  KKT(P)  |  KKT(D)  | KKT(D,s) |  KKT(C)  | KKT(C,s) |  KKT(T)  | VI_nat_res |  time(s) |';
    disp(headMsg)
    stage_1_Msg = ['  ',...
        num2str(0,'%10.2d'), '/', num2str(continuationStepNum,'%10.2d'),'  | ',...
        num2str(Log.param(1, 1), '%10.2e'),' | ',...
        num2str(Log.param(1, 2), '%10.2e'),' | ',...
        num2str(Log.cost(1),'%10.2e'), ' | ',...
        num2str(Log.KKT_error(1, 1), '%10.2e'), ' | ',...
        num2str(Log.KKT_error(1, 2), '%10.2e'), ' | ',...
        num2str(Log.KKT_error(1, 3), '%10.2e'), ' | ',...
        num2str(Log.KKT_error(1, 4), '%10.2e'), ' | ',...
        num2str(Log.KKT_error(1, 5), '%10.2e'), ' | ',...
        num2str(Log.KKT_error(1, 6), '%10.2e'), ' |  ',...
        num2str(Log.VI_natural_residual(1), '%10.2e'),'  |  ',...
        num2str(Log.timeElapsed(1), '%10.4f'),'  | '];
    disp(stage_1_Msg)    
    
    %% check ternimation based on stage 1's terminal status 
    if Info_Stage_1.terminalStatus ~= 1
        % stage 1 fails: terminal the overall algorithm
        % -- failure case 1: non-interior-point method fails to solve the first NLP
        terminalStatus = 0;
        terminalMsg = Info_Stage_1.terminalMsg;
        break
    else
        % stage 1 success: check whether need to performs Euler-Newton continuation method
        if continuationStepNum == 0
            % -- success case 1: the first NLP is solved and the given parameter does not allow Euler-Newton continuation method
            terminalStatus = 1;
            terminalMsg = ['- Solver succeeds: ',...
                'because the first NLP is solved by NIP method and the given parameter does not allow the continuation method'];
            break
        end
    end

    %% stage 2: solve the subsequent parameterized NLPs by Euler-Newton continuation method
    % estimate a new parameter for the first Euler-Newton continuation step
    j = 1;
    p_trial = min([kappa_times .* p, p.^kappa_exp], [], 2);
    p_j = max([p_trial, p_End], [], 2);
    % iteration routine (z: previous iterate z_{j-1}, z_j: current iterate z_{j})
    while true
        %% Euler predict step
        timeStart_continuationStep = tic;
        
        [gamma_h_Eu_j, gamma_c_Eu_j, gamma_g_Eu_j, z_Eu_j, Info_Euler] =...
            self.Euler_predict_step(gamma_h, gamma_c, gamma_g, z, p, p_j);
        Time.Euler_predict_funcGradEval = Time.Euler_predict_funcGradEval + Info_Euler.timeElasped_funcGradEval;
        Time.Euler_predict_stepEval = Time.Euler_predict_stepEval + Info_Euler.timeElasped_stepEval;
        
        %% Newton correct step
        [gamma_h_j, gamma_c_j, gamma_g_j, z_j, Info_Newton] = ...
            self.Newton_correct_step(gamma_h_Eu_j, gamma_c_Eu_j, gamma_g_Eu_j, z_Eu_j, p_j);
        Time.Newton_correct_funcGradEval = Time.Newton_correct_funcGradEval + Info_Newton.timeElasped_funcGradEval;
        Time.Newton_correct_stepEval = Time.Newton_correct_stepEval + Info_Newton.timeElasped_stepEval;
        
        timeElapsed_continuationStep = toc(timeStart_continuationStep);
        Time.total = Time.total + timeElapsed_continuationStep;
        if Info_Newton.terminalStatus ~= 1
            EulerNewton_terminalStatus = 0; % Newton step fails
            break
        end
        
        %% evaluate cost, KKT error, and VI natural residual
        % cost, constraint and their Jacobian
        J_j = full(NLP.FuncObj.J(z_j));
        h_j = full(NLP.FuncObj.h(z_j));
        c_j = full(NLP.FuncObj.c(z_j));
        g_j = full(NLP.FuncObj.g(z_j, p_j(1)));
        J_grad_j = full(NLP.FuncObj.J_grad(z_j));
        h_grad_j = sparse(NLP.FuncObj.h_grad(z_j));
        c_grad_j = sparse(NLP.FuncObj.c_grad(z_j));
        g_grad_j = sparse(NLP.FuncObj.g_grad(z_j, p_j(1)));
        % KKT error
        LAG_grad_z_j = J_grad_j + gamma_h_j' * h_grad_j - gamma_c_j' * c_grad_j  - gamma_g_j' * g_grad_j;
        scaling_dual_j = max([Option.KKT_scaling_max,...
            norm([gamma_h_j; gamma_c_j; gamma_g_j], 1)/(NLP.Dim.h + NLP.Dim.c + NLP.Dim.g)])/Option.KKT_scaling_max;
        scaling_complementary_j = max([Option.KKT_scaling_max,...
            norm([gamma_c_j; gamma_g_j], 1)/(NLP.Dim.c + NLP.Dim.g)])/Option.KKT_scaling_max;
        KKT_error_primal_j = norm([h_j;...
            min([zeros(NLP.Dim.c, 1), c_j], [], 2);...
            min([zeros(NLP.Dim.g, 1), g_j], [], 2)], inf);
        KKT_error_dual_j = norm([LAG_grad_z_j';...
            min([zeros(NLP.Dim.c, 1), gamma_c_j], [], 2);...
            min([zeros(NLP.Dim.g, 1), gamma_g_j], [], 2)], inf);
        KKT_error_dual_scaled_j = KKT_error_dual_j/scaling_dual_j;
        KKT_error_complementary_j = norm([c_j .* gamma_c_j; g_j .* gamma_g_j], inf);
        KKT_error_complementary_scaled_j = KKT_error_complementary_j/scaling_complementary_j;
        KKT_error_total_j = max([KKT_error_primal_j, KKT_error_dual_scaled_j, KKT_error_complementary_scaled_j]);
        % VI natural residual
        VI_natural_residual_j  = self.evaluateNaturalResidual(z_j);
        
        %% record and print information of this iteration j
        % record
        Log.param(j + 1, :) = p_j';
        Log.cost(j + 1, :)  = J_j;
        Log.KKT_error(j + 1, :) = [KKT_error_primal_j, KKT_error_dual_j, KKT_error_dual_scaled_j,...
            KKT_error_complementary_j, KKT_error_complementary_scaled_j, KKT_error_total_j];
        Log.VI_natural_residual(j + 1, :) = VI_natural_residual_j;
        Log.timeElapsed(j + 1, :) = timeElapsed_continuationStep;
        % print
        if mod(j, 10) ==  0
            disp('---------------------------------------------------------------------------------------------------------------------------------------------')
            headMsg = ' StepNum |    s     |   sigma  |   cost   |  KKT(P)  |  KKT(D)  | KKT(D,s) |  KKT(C)  | KKT(C,s) |  KKT(T)  | VI_nat_res |  time(s) |';
            disp(headMsg)
        end        
        continuation_Step_Msg = ['  ',...
            num2str(j,'%10.2d'), '/', num2str(continuationStepNum,'%10.2d'),'  | ',...
            num2str(Log.param(j + 1, 1), '%10.2e'),' | ',...
            num2str(Log.param(j + 1, 2), '%10.2e'),' | ',...
            num2str(Log.cost(j + 1),'%10.2e'), ' | ',...
            num2str(Log.KKT_error(j + 1, 1), '%10.2e'), ' | ',...
            num2str(Log.KKT_error(j + 1, 2), '%10.2e'), ' | ',...
            num2str(Log.KKT_error(j + 1, 3), '%10.2e'), ' | ',...
            num2str(Log.KKT_error(j + 1, 4), '%10.2e'), ' | ',...
            num2str(Log.KKT_error(j + 1, 5), '%10.2e'), ' | ',...
            num2str(Log.KKT_error(j + 1, 6), '%10.2e'), ' |  ',...
            num2str(Log.VI_natural_residual(j + 1), '%10.2e'),'  |  ',...
            num2str(Log.timeElapsed(j + 1), '%10.4f'),'  | '];
        disp(continuation_Step_Msg)
        
        %% prepare next iteration
        if j == continuationStepNum
            % Euler Newton continuation step has all been done, break
            EulerNewton_terminalStatus = 1;
            break
        else
            j = j + 1;
            % iterate and parameter
            gamma_h = gamma_h_j;
            gamma_c = gamma_c_j;
            gamma_g = gamma_g_j;
            z = z_j;
            p = p_j;
            p_trial = min([kappa_times .* p, p.^kappa_exp], [], 2);
            p_j = max([p_trial, p_End], [], 2);
            % cost, KKT error, and VI residual
            J = J_j;
            KKT_error_primal               = KKT_error_primal_j;
            KKT_error_dual                 = KKT_error_dual_j;
            KKT_error_dual_scaled          = KKT_error_dual_scaled_j;
            KKT_error_complementary        = KKT_error_complementary_j;
            KKT_error_complementary_scaled = KKT_error_complementary_scaled_j;
            KKT_error_total                = KKT_error_total_j;
            VI_natural_residual            = VI_natural_residual_j;
        end
        
    end
    
    %% check termination based on Euler-Newton continuation method' terminal status
    disp('Done!')
    disp('---------------------------------------------------------------------------------------------------------------------------------------------')
    switch EulerNewton_terminalStatus
        case 0
            terminalStatus = 0;
            terminalMsg = ['- Solver fails: ',...
                'because Euler-Newton continuation method fails to find a stepsize'];
        case 1
            terminalStatus = 1;
            terminalMsg = ['- Solver succeeds: ',...
                'because Euler-Newton continuation method has been done'];
    end
    break
end

%% return optimal solution and create information
% return previous iterate as solution
z_Opt = z;
% create Info (basic: time, continuation step  terminal status)
timeElapsed_Euler_Newton = Time.Euler_predict_funcGradEval + Time.Euler_predict_stepEval...
    + Time.Newton_correct_funcGradEval + Time.Newton_correct_stepEval;
Time.else = Time.total - Time.non_interior_point - timeElapsed_Euler_Newton;
Info.Time = Time;
Info.continuationStepNum = continuationStepNum;
Info.terminalStatus = terminalStatus;
Info.terminalMsg = terminalMsg;
% create Info (stage 1)
Info.Stage_1 = Info_Stage_1;
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
Info.VI_natural_residual            = VI_natural_residual;
% create Info (log)
Info.Log = Log;
% display termination and solution message
disp('*--------------------------------------------- Solution Information ----------------------------------------------*')
disp('1. Terminal Status')
disp(Info.terminalMsg)
disp('2. Continuation Step Message')
disp(['- Continuation Step: ............................', num2str(Info.continuationStepNum)])
disp(['- TimeElapsed: ................................. ', num2str(Info.Time.total,'%10.3f'), 's'])
disp(['- Time for Stage 1: .............................', num2str(Info.Time.non_interior_point,'%10.3f'), 's'])
disp(['- Time for Euler Newton Step: ...................', num2str(timeElapsed_Euler_Newton,'%10.3f'), 's'])
disp(['- Average Time for Euler Newton Step: ...........', num2str(1000 * (timeElapsed_Euler_Newton)/Info.continuationStepNum,'%10.2f'), ' ms/step'])
disp('3. Solution Message')
disp(['- Cost: ........................................ ', num2str(Info.cost,'%10.3e'), '; '])
disp(['- KKT (primal): ................................ ', num2str(Info.KKT_error_primal,'%10.3e'), '; '])
disp(['- KKT (dual, unscaled / scaled): ............... ', num2str(Info.KKT_error_dual,'%10.3e'), ' / ', num2str(Info.KKT_error_dual_scaled,'%10.3e')  '; '])
disp(['- KKT (complementary, unscaled / scaled): ...... ', num2str(Info.KKT_error_complementary,'%10.3e'), ' / ', num2str(Info.KKT_error_complementary_scaled,'%10.3e'), '; '])
disp(['- KKT (total): ................................. ', num2str(Info.KKT_error_total,'%10.3e'), '; '])
disp(['- VI natural residual: ......................... ', num2str(Info.VI_natural_residual,'%10.3e'), '; '])
end

