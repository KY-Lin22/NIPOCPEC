function [z_Opt, Info] = solveNLP_IP_homotopy(OCPEC, NLP, solver, z_Init, s_Init, s_End)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*

% check input z_Init
if ~all(size(z_Init) == [NLP.Dim.z, 1])
    error('z_Init has wrong dimension')
end

% check scalar negative parameter
if (~all(size(s_Init) == [1, 1])) || (s_Init < 0)
    error('s_Init should be scalar and nonnegative')
end
if (~all(size(s_End) == [1, 1])) || (s_End < 0)
    error('s_End should be scalar and nonnegative')
end
if s_Init < s_End
    error('s_Init should not smaller than s_End')
end

% homotopy parameter
kappa_s_times = 0.9;
kappa_s_exp = 1.1;

% evaluate max homotopy counter based on s_Init and s_End
s_test = s_Init;
maxHomotopyIter = 1;
while true
    if s_test == s_End
        break
    else
        s_trial = min([kappa_s_times .* s_test, s_test.^kappa_s_exp]);
        s_test = max([s_trial, s_End]);
        maxHomotopyIter = maxHomotopyIter + 1;
    end
end

% log (time, param, cost, KKT error, natRes)
Log.param      = zeros(maxHomotopyIter, 1); % s
Log.cost       = zeros(maxHomotopyIter, 1);
Log.KKT_error  = zeros(maxHomotopyIter, 3); % [primal, dual_scaled, total]
Log.VI_nat_res = zeros(maxHomotopyIter, 1);
Log.iterNum    = zeros(maxHomotopyIter, 1);
Log.time       = zeros(maxHomotopyIter, 1); % elapsed time in each continuation step

%% homotopy loop (j: homotopy counter)
z_Init_j = z_Init;
s_j = s_Init;
timeElasped = 0;
for j = 1 : maxHomotopyIter
    % step 1: solve a NLP with given s
    timeStart_j = tic;
    solution_j = solver('x0', z_Init_j, 'p', s_j,...
        'lbg', [zeros(NLP.Dim.h, 1); zeros(NLP.Dim.c, 1); zeros(NLP.Dim.g, 1)],...
        'ubg', [zeros(NLP.Dim.h, 1); inf*ones(NLP.Dim.c, 1); inf*ones(NLP.Dim.g, 1)]);
    z_Opt_j = full(solution_j.x);
    VI_nat_res_j = norm(evaluateNaturalResidual_IP(OCPEC, NLP, z_Opt_j), inf);
    time_j = toc(timeStart_j); 
    
    % step 2: record and print information of the current homotopy iterate
    KKT_error_primal_j = solver.stats.iterations.inf_pr(end);
    KKT_error_dual_j = solver.stats.iterations.inf_du(end);
    timeElasped = timeElasped + time_j;
    Log.param(j) = s_j;
    Log.cost(j) = full(solution_j.f);
    Log.KKT_error(j, :) = [KKT_error_primal_j, KKT_error_dual_j, max(KKT_error_primal_j, KKT_error_dual_j)];
    Log.VI_nat_res(j) = VI_nat_res_j;
    Log.iterNum(j) = solver.stats.iter_count;
    Log.time(j) = time_j;
    
    if mod(j, 10) == 1
        disp('---------------------------------------------------------------------------------------------------------------------------------------------')
        headMsg = ' Homotopy |    s     |   cost   |  KKT(P)  |  KKT(D)  | VI_nat_res | iterNum |  time(s) |';
        disp(headMsg)
    end
    prevIterMsg = ['  ',...
        num2str(j,'%10.2d'), '/', num2str(maxHomotopyIter,'%10.2d'),'   | ',...
        num2str(Log.param(j), '%10.2e'),' | ',...
        num2str(Log.cost(j), '%10.2e'),' | ',...
        num2str(Log.KKT_error(j, 1), '%10.2e'),' | ',...
        num2str(Log.KKT_error(j, 2), '%10.2e'),' |  ',...
        num2str(Log.VI_nat_res(j), '%10.2e'),'  |   ',...
        num2str(Log.iterNum(j), '%10.3d'),'   |  ',...
        num2str(Log.time(j), '%10.4f'),'  | '];
    disp(prevIterMsg)

    % step 3: check ternimation based on the current homotopy iterate
    if strcmp(solver.stats.return_status, 'Solve_Succeeded') && (j == maxHomotopyIter)
        % IPOPT at the final homotopy iteration finds the optimal solution
        exitFlag = true;
    elseif ~strcmp(solver.stats.return_status, 'Solve_Succeeded')
        % IPOPT at this homotopy iteration fails to find the optimal solution
        exitFlag = true;
    else
        % IPOPT at this homotopy iteration (not the final) finds the optimal solution, prepare for next homotopy iteration
        exitFlag = false;
        z_Init_j = z_Opt_j;
        s_trial = min([kappa_s_times .* s_j, s_j.^kappa_s_exp]);
        s_j = max([s_trial, s_End]);        
    end
    
    % step 4: check exitFlag and return optimal solution
    if exitFlag
       % return the current homotopy iterate as the optimal solution
       z_Opt = z_Opt_j; 
       % create Info
       Info.time = timeElasped;
       Info.cost = full(solution_j.f);
       Info.VI_natural_residual = VI_nat_res_j;
       Info.Log = Log;
       % display homotopy terminal result and then break
       disp('Done!')
       disp(solver.stats.return_status)
       disp(['Time: ', num2str(Info.time,'%10.3f'), ' s; ',...
           'Cost: ', num2str(Info.cost,'%10.5f'), '; ',...
           'VI natural residual: ', num2str(Info.VI_natural_residual,'%10.3e')])
       
       break
    end
end

end
