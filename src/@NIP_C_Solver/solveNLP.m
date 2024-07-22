function [z_Opt, Info] = solveNLP(self, z_Init)
% solve NLP with given z_Init by non-interior-point continuation method
% NLP has the form:
%  min  J(z),
%  s.t. h(z) = 0,
%       c(z) >= 0,
%       g(z, s) >= 0
% where: z is the variable,
%        s is the parameter,
%        J is the cost, and h, c, g are the constraint
% Syntax:
%          [z_Opt, Info] = solveNLP(self, z_Init)
%          [z_Opt, Info] = self.solveNLP(z_Init)
% Argument:
%          z_Init: double, NLP.Dim.z X 1, initial guess
% Output:
%          z_Opt: double, NLP.Dim.z X 1, optimal solution found by solver
%          Info: struct, record the iteration information
import casadi.*

%% check input
if ~all(size(z_Init) == [self.NLP.Dim.z, 1])
    error('z_Init has wrong dimension')
end

%% Initialization
% Y node (gamma_h, gamma_c, gamma_g, z)
Y_Node = cumsum([self.NLP.Dim.h, self.NLP.Dim.c, self.NLP.Dim.g, self.NLP.Dim.z]);
% preproccess initial guess
z_Init = self.preprocessInitialGuess(z_Init);
% generate parameter sequence
[P, l_Max] = self.createParameterSequence();
% create record
Log.param          = zeros(l_Max + 1, 2); % [s, sigma]
Log.cost       = zeros(l_Max + 1, 1);
Log.KKT_res    = zeros(l_Max + 1, 1);
Log.KKT_error  = zeros(l_Max + 1, 1); 
Log.VI_nat_res = zeros(l_Max + 1, 1);
Log.time       = zeros(l_Max + 1, 1); 

%% continuation loop (l: continuation step counter, Y_l: current iterate, Y: previous iterate)
l = 0;
% continuously update p_l and Y_l
while true
    %% step 1: evaluate iterate at current continuation step
    % specify parameter p
    p_l = P(:, l + 1);
    % evaluate iterate Y
    if l == 0
        % stage 1: initialize iterate Y_l by non-interior-point method
        [Y_l, Info_NIP] = self.NonInteriorPointMethod(z_Init, p_l);
        terminal_status_l = Info_NIP.terminal_status;
        terminal_msg_l = Info_NIP.terminal_msg;
        timeElasped_l = Info_NIP.time;        
    else
        % stage 2: evaluate iterate Y_l by Euler-Newton (predictor-corrector) continuation method
        [Y_l, Info_EulerNewton] = self.EulerNewtonContinuationMethod(Y, p_l, p);
        terminal_status_l = Info_EulerNewton.terminal_status;
        terminal_msg_l = Info_EulerNewton.terminal_msg; 
        timeElasped_l = Info_EulerNewton.time;        
    end

    %% step 2: record and print information of the current continuation step
    % KKT residual (FB reformulation)
    KKT_residual_l = full(self.FuncObj.KKT_residual(Y_l, p_l)); 
    % some NLP quantities
    z_l = Y_l(Y_Node(3) + 1 : Y_Node(4), 1);    
    J_l = full(self.FuncObj.J(z_l));
    [KKT_error_primal_l, KKT_error_dual_l, KKT_error_complementary_l, KKT_error_total_l] = self.evaluateKKT_Error(Y_l, p_l);
    VI_nat_res_l = norm(self.evaluateNaturalResidual(z_l), inf);
    % record
    Log.param(l + 1, :) = p_l';  
    Log.cost(l + 1, :) = J_l;
    Log.KKT_res(l + 1, :) = norm(KKT_residual_l, 2)/self.OCPEC.nStages;
    Log.KKT_error(l + 1, :) = KKT_error_total_l;
    Log.VI_nat_res(l + 1, :) = VI_nat_res_l;
    Log.time(l + 1, :) = timeElasped_l;
    % print
    if mod(l, 10) ==  0
        disp('--------------------------------------------------------------------------------------------------------------------------------')
        headMsg = '   StepNum  |  param(s / sigma) |   cost   | KKT_res  | KKT_error | KKT_error(Pri / Dual / Comp) | VI_nat_res | time[s] ';
        disp(headMsg)
    end
    continuation_Step_Msg = ['  ',...
        num2str(l,'%10.4d'), '/', num2str(l_Max,'%10.4d'),' | ',...
        num2str(Log.param(l + 1, 1), '%10.2e'),' ', num2str(Log.param(l + 1, 2), '%10.2e'),' | ',...
        num2str(Log.cost(l + 1),'%10.2e'), ' | ',...
        num2str(Log.KKT_res(l + 1), '%10.2e'), ' |  ',...
        num2str(Log.KKT_error(l + 1), '%10.2e'), ' |  ',...
        num2str(KKT_error_primal_l, '%10.2e'), ' ', num2str(KKT_error_dual_l, '%10.2e'), ' ', num2str(KKT_error_complementary_l, '%10.2e'), '  |  ',...
        num2str(Log.VI_nat_res(l + 1), '%10.2e'),'  | ',...
        num2str(Log.time(l + 1), '%10.4f')];
    disp(continuation_Step_Msg)

    %% step 3: check ternimation based on the current continuation step
    if terminal_status_l && (VI_nat_res_l <= self.Option.Continuation.tol.VI_nat_res) && (KKT_error_l <= self.Option.Continuation.tol.KKT_error)
        % this continuation step finds the desired optimal solution
        terminal_status = 1;
        terminal_msg = terminal_msg_l;
        break
    elseif ~terminal_status_l
        % this continuation step fails to find the optimal solution
        terminal_status = 0;
        terminal_msg = terminal_msg_l;
        break
    elseif l == l_Max
        % final continuation step still can not find the desired optimal solution
        terminal_status = -1;
        terminal_msg = 'solver can not find the optimal solution satisfying the desired VI natural residual';
        break
    else
        % this continuation step finds the optimal solution, prepare for next step
        l = l + 1;
        Y = Y_l;
        p = p_l;
    end

end

%% return optimal solution and create information
% extract primal and dual variable
gamma_h_l = Y_l(            1 : Y_Node(1), 1);
gamma_c_l = Y_l(Y_Node(1) + 1 : Y_Node(2), 1);
gamma_g_l = Y_l(Y_Node(2) + 1 : Y_Node(3), 1);
z_l       = Y_l(Y_Node(3) + 1 : Y_Node(4), 1);
% return the current homotopy iterate as the optimal solution
z_Opt = z_l;
% create Info
Info.continuationStepNum = l;
Info.terminal_status = terminal_status;
Info.terminal_msg = terminal_msg;
Info.gamma_h = gamma_h_l;
Info.gamma_c = gamma_c_l;
Info.gamma_g = gamma_g_l;
Info.cost = J_l;
Info.KKT_error = KKT_error_total_l;
Info.VI_natural_residual = VI_nat_res_l;
Info.Log = Log;
Info.time = sum(Log.time);
% display result
disp('*------------------------------------------------- Solution Information --------------------------------------------------*')
disp('1. Terminal Message')
disp(Info.terminal_msg)
disp('2. Continuation Step Message')
disp(['- TimeElapsed: ................................................... ', num2str(Info.time,'%10.4f'), ' s'])
disp(['- Total time for solving first parameterized NLP: ................ ', num2str(Info.Log.time(1),'%10.4f'), ' s'])
disp(['- Total time for solving subsequent parameterized NLP: ........... ', num2str((Info.time - Info.Log.time(1)),'%10.4f'), ' s'])
disp(['- Continuation Step: ............................................. ', num2str(Info.continuationStepNum)])
if Info.continuationStepNum ~= 0
    disp(['- Average time for each step using Euler Newton method: ...... ', num2str((Info.time - Info.Log.time(1))/Info.continuationStepNum,'%10.4f'), ' s'])
end
disp('3. Solution Message')
disp(['- Cost: .......................................................... ', num2str(Info.cost,'%10.3e'), '; '])
disp(['- KKT error: ..................................................... ', num2str(Info.KKT_error,'%10.3e'), '; '])
disp(['- VI natural residual: ........................................... ', num2str(Info.VI_natural_residual,'%10.3e'), '; '])


end