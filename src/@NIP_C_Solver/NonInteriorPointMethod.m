function [Y_Opt, Info] = NonInteriorPointMethod(self, z_Init, p_Init)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

%% iteration routine (Y: previous iterate Y_{k-1}, Y_k: current iterate Y_{k}) 
% time record
Time = struct('searchDirection', 0, 'lineSearch', 0, 'else', 0, 'total', 0);
% counter, beta and iterate
k = 1;
beta = self.Option.LineSearch.betaInit;
Y_Node = cumsum([self.NLP.Dim.h, self.NLP.Dim.c, self.NLP.Dim.g, self.NLP.Dim.z]);
Y = [ones(self.NLP.Dim.h, 1); ones(self.NLP.Dim.c, 1); ones(self.NLP.Dim.g, 1); z_Init];
% iteration loop
while true
    %% step 0: check iteration counter
    if k > self.Option.maxIterNum
        % failure case 1: exceed the maximum number of iteration
        terminal_status = 0;
        terminal_msg = ['- Solver fails: ', 'because the maximum number of iteration exceeded'];
        break
    end
    timeStart_total = tic;

    %% step 1: search direction evaluation based on previous iterate z
    timeStart_SearchDirection = tic;    
    dY = full(self.FuncObj.dY_Newton(Y, p_Init));         
    timeElasped_searchDirection = toc(timeStart_SearchDirection);
    
    %% step 2: check whether we can terminate successfully based on the previous iterate
    [KKT_error_primal, KKT_error_dual, KKT_error_complementary, KKT_error_total] = self.evaluateKKT_Error(Y, p_Init);
    dYNorm = norm(dY, inf);% dYNorm (L_inf norm) 
    if KKT_error_total < self.Option.tol.KKT_error_total
        % Success case 1: the KKT error satisfies tolerance
        terminal_status = 1;
        terminal_msg = ['- Solver succeeds: ', 'because the KKT error satisfies tolerance'];
        break
    elseif dYNorm < self.Option.tol.dYNorm
        % Success case 2: the norm of search direction satisfies tolerance
        terminal_status = 1;
        terminal_msg = ['- Solver succeeds: ', 'because the norm of search direction satisfies tolerance'];
        break
    elseif (KKT_error_primal <= self.Option.tol.KKT_error_primal) && ...
            (KKT_error_dual <= self.Option.tol.KKT_error_dual) && ...
            (KKT_error_complementary <= self.Option.tol.KKT_error_complementarity)
        % Success case 3: primal and dual error satisfy individual tolerance
        terminal_status = 1;
        terminal_msg = ['- Solver succeeds: ', 'because primal, dual, and complementarity error satisfy individual tolerance'];
        break
    end
    %% step 3: merit line search
    [Y_k, Info_LineSearch] = self.LineSearchMerit(beta, Y, p_Init, dY);
    % check status
    if Info_LineSearch.status == 0
        % failure case 2: line search fails
        terminal_status = 0;
        terminal_msg = ['- Solver fails: ', 'because merit line search reaches the min stepsize'];        
        break
    else
        % line search quantities
        beta_k   = Info_LineSearch.beta;
        stepSize = Info_LineSearch.stepSize;
        merit    = Info_LineSearch.merit;
    end
    timeElasped_lineSearch = Info_LineSearch.time;

    %% step 4: record and print information of this iteration k
    timeElasped_total = toc(timeStart_total);
    Time.searchDirection = Time.searchDirection + timeElasped_searchDirection;
    Time.lineSearch = Time.lineSearch + timeElasped_lineSearch;
    Time.total = Time.total + timeElasped_total;
    % some other quantities
    z = Y(Y_Node(3) + 1 : Y_Node(4), 1);
    J = full(self.FuncObj.J(z));
    % head
    if mod(k, 10) == 1
        disp('----------------------------------------------------------------------------------------------------------------------')
        headMsg = ' Iter |   cost   |  KKT(P)  |  KKT(D)  |  KKT(C)  |  dYNorm  |   beta   | stepsize |  merit   | merit(t) | time(ms) |';
        disp(headMsg)
    end
    % previous iterate message
    prevIterMsg = ['  ',...
        num2str(k,'%10.3d'),' | ',...
        num2str(J,'%10.2e'), ' | ',...
        num2str(KKT_error_primal, '%10.2e'), ' | ',...
        num2str(KKT_error_dual, '%10.2e'),  ' | ',...
        num2str(KKT_error_complementary, '%10.2e'), ' | ',...
        num2str(dYNorm,'%10.2e'), ' | ',...
        num2str(beta_k,'%10.2e'), ' | ',...
        num2str(stepSize,'%10.2e'), ' | ', ...
        num2str(merit(1),'%10.2e'), ' | ', num2str(merit(2),'%10.2e'), ' | ',...
        num2str(1000 * timeElasped_total,'%10.2e'), ' | '];
    disp(prevIterMsg)

    %% step 6: prepare next iteration
    k = k + 1;
    beta = beta_k;
    Y = Y_k;    

end

%% return optimal solution and create information
% return previous iterate as solution
Y_Opt = Y;
z = Y(Y_Node(3) + 1 : Y_Node(4), 1);
% create Info (basic: time, iterNum, terminalStatus)
Time.else = Time.total - Time.searchDirection - Time.lineSearch;
Info.Time = Time;
Info.time = Time.total;
Info.iterNum = k - 1;
Info.terminal_status = terminal_status;
Info.terminal_msg = terminal_msg;
% create Info (corresponds to the solution: cost, KKT, natural residual)
Info.cost = J;
Info.KKT_error = KKT_error_total;
Info.VI_natural_residual = norm(self.evaluateNaturalResidual(z), inf);
% display termination and solution message
disp('*--------------------------------------------- Solution Information ----------------------------------------------*')
disp('1. Terminal Status')
disp(Info.terminal_msg)
disp('2. Iteration Process Message')
disp(['- Iterations: ................... ', num2str(Info.iterNum)])
disp(['- TimeElapsed: .................. ', num2str(Info.Time.total,'%10.3f'), 's'])
disp(['- AverageTime: .................. ', num2str(1000 * Info.Time.total /Info.iterNum,'%10.2f'), ' ms/Iter'])
disp('3. Solution Message')
disp(['- Cost: ......................... ', num2str(Info.cost,'%10.3e'), '; '])
disp(['- KKT error: .................... ', num2str(Info.KKT_error, '%10.3e'), '; '])
disp(['- VI natural residual: .......... ', num2str(Info.VI_natural_residual,'%10.3e'), '; '])
    
end