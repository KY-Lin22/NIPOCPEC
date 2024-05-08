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
%% check input
% check input z_Init
if ~all(size(z_Init) == [self.NLP.Dim.z, 1])
    error('z_Init has wrong dimension')
end
% check relaxation parameter
if (p_Init(1) < 0) || (p_End(1) < 0)
    error('relax parameter s (i.e., p_1) should be nonnegative')
end
if p_Init(1) < p_End(1)
    error('s_Init should not smaller than s_End')
end
% check FB smooth parameter
if (p_Init(2) < 0) || (p_End(2) < 0)
    error('FB smooth parameter sigma (i.e, p_2) should be nonnegative')
end
if p_Init(2) < p_End(2)
    error('sigma_Init should not smaller than sigma_End')
end
% load parameter
kappa_times = self.Option.Homotopy.kappa_times;
kappa_exp = self.Option.Homotopy.kappa_exp;
% preproccess initial guess
z_Init = self.preprocessInitialGuess(z_Init);

%% create record for log 
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
% log
Log.param               = zeros(continuationStepNum + 1, 2); % [s, sigma]
Log.cost                = zeros(continuationStepNum + 1, 1);
Log.KKT_error           = zeros(continuationStepNum + 1, 6); % [primal, dual, dual_scaled, complementary, complementary_scaled, total]
Log.VI_natural_residual = zeros(continuationStepNum + 1, 1);
Log.timeElapsed         = zeros(continuationStepNum + 1, 1); % elapsed time in each continuation step
% time
Time_continuation = struct('Euler_funcEval', 0, 'Euler_stepEval', 0,'Newton_funcEval', 0, 'Newton_stepEval',0, 'total', 0);

%% non-interior-point continuation method (z: previous iterate z_{j-1}, z_j: current iterate z_{j})
while true
    %% stage 1: solve the first parameterized NLP by non-interior-point method (compute z based on p)
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
    Time_non_interior_point        = Info_Stage_1.Time;
    % record          
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
    % check ternimation based on stage 1's terminal status 
    if Info_Stage_1.terminalStatus ~= 1
        % stage 1 fails: non-interior-point method fails to solve the first NLP, terminal the overall algorithm
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
    for j = 1 : continuationStepNum
        % estimate a new parameter for the Euler-Newton continuation step
        p_trial = min([kappa_times .* p, p.^kappa_exp], [], 2);
        p_j = max([p_trial, p_End], [], 2);
        % Euler predict step
        [gamma_h_Eu_j, gamma_c_Eu_j, gamma_g_Eu_j, z_Eu_j, Info_Euler] = self.Euler_predict_step(gamma_h, gamma_c, gamma_g, z, p, p_j);
        % Newton correct step
        [gamma_h_j, gamma_c_j, gamma_g_j, z_j, Info_Newton] = self.Newton_correct_step(gamma_h_Eu_j, gamma_c_Eu_j, gamma_g_Eu_j, z_Eu_j, p_j);
        % record time
        Time_continuation.Euler_funcEval  = Time_continuation.Euler_funcEval  + Info_Euler.time.funcEval;
        Time_continuation.Euler_stepEval  = Time_continuation.Euler_stepEval  + Info_Euler.time.stepEval;
        Time_continuation.Newton_funcEval = Time_continuation.Newton_funcEval + Info_Newton.time.funcEval;
        Time_continuation.Newton_stepEval = Time_continuation.Newton_stepEval + Info_Newton.time.stepEval;
        timeElapsed_j = Info_Euler.time.funcEval + Info_Euler.time.stepEval + Info_Newton.time.funcEval + Info_Newton.time.stepEval;
        Time_continuation.total = Time_continuation.total + timeElapsed_j;
        % extract quantities (cost, KKT error, VI residual)
        J_j = Info_Newton.cost;
        KKT_error_primal_j               = Info_Newton.KKT_error_primal;
        KKT_error_dual_j                 = Info_Newton.KKT_error_dual;
        KKT_error_dual_scaled_j          = Info_Newton.KKT_error_dual_scaled;
        KKT_error_complementary_j        = Info_Newton.KKT_error_complementary;
        KKT_error_complementary_scaled_j = Info_Newton.KKT_error_complementary_scaled;
        KKT_error_total_j                = Info_Newton.KKT_error_total;
        VI_natural_residual_j            = Info_Newton.VI_natural_residual;
        % record
        Log.param(j + 1, :) = p_j';
        Log.cost(j + 1, :)  = J_j;
        Log.KKT_error(j + 1, :) = [KKT_error_primal_j, KKT_error_dual_j, KKT_error_dual_scaled_j,...
            KKT_error_complementary_j, KKT_error_complementary_scaled_j, KKT_error_total_j];
        Log.VI_natural_residual(j + 1, :) = VI_natural_residual_j;
        Log.timeElapsed(j + 1, :) = timeElapsed_j;
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
        % prepare next iteration: iterate and parameter
        gamma_h = gamma_h_j;
        gamma_c = gamma_c_j;
        gamma_g = gamma_g_j;
        z = z_j;
        p = p_j;
        % prepare next iteration: cost, KKT error, and VI residual
        J = J_j;
        KKT_error_primal               = KKT_error_primal_j;
        KKT_error_dual                 = KKT_error_dual_j;
        KKT_error_dual_scaled          = KKT_error_dual_scaled_j;
        KKT_error_complementary        = KKT_error_complementary_j;
        KKT_error_complementary_scaled = KKT_error_complementary_scaled_j;
        KKT_error_total                = KKT_error_total_j;
        VI_natural_residual            = VI_natural_residual_j;
    end  
    % terminate Euler-Newton continuation method
    terminalStatus = 1;
    terminalMsg = ['- Solver succeeds: ', 'because Euler-Newton continuation method has been done'];
    break
end

%% return optimal solution and create information
% return previous iterate as solution
z_Opt = z;
% create Info (basic: time, continuation step  terminal status)
Time.non_interior_point = Time_non_interior_point;
Time.continuation = Time_continuation;
Time.total = Time_non_interior_point.total + Time_continuation.total;
Info.Time = Time;
Info.continuationStepNum = continuationStepNum;
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
disp(['- Time for Stage 1: .............................', num2str(Info.Time.non_interior_point.total,'%10.3f'), 's'])
disp(['- Time for Euler Newton Step: ...................', num2str(Info.Time.continuation.total,'%10.3f'), 's'])
disp(['- Average Time for Euler Newton Step: ...........', num2str(1000 * (Info.Time.continuation.total)/Info.continuationStepNum,'%10.2f'), ' ms/step'])
disp('3. Solution Message')
disp(['- Cost: ........................................ ', num2str(Info.cost,'%10.3e'), '; '])
disp(['- KKT (primal): ................................ ', num2str(Info.KKT_error_primal,'%10.3e'), '; '])
disp(['- KKT (dual, unscaled / scaled): ............... ', num2str(Info.KKT_error_dual,'%10.3e'), ' / ', num2str(Info.KKT_error_dual_scaled,'%10.3e')  '; '])
disp(['- KKT (complementary, unscaled / scaled): ...... ', num2str(Info.KKT_error_complementary,'%10.3e'), ' / ', num2str(Info.KKT_error_complementary_scaled,'%10.3e'), '; '])
disp(['- KKT (total): ................................. ', num2str(Info.KKT_error_total,'%10.3e'), '; '])
disp(['- VI natural residual: ......................... ', num2str(Info.VI_natural_residual,'%10.3e'), '; '])
end

