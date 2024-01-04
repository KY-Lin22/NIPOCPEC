function showResult(self, Info)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Option = self.Option;
%% Stage 1: Log and time
% plot iteration log
if Option.recordLevel == 1
    stage_1_Log = Info.Stage_1.Log;
    stage_1_logAxis = 0 : 1 : size(stage_1_Log.dYNorm, 1) - 1; % define log axis       
    figure(101) 
    % dY norm
    subplot(4,1,1)
    plot(stage_1_logAxis, log10(stage_1_Log.dYNorm))
    ylabel('log10')
    title('stage 1: norm of search direction')
    % step size
    subplot(4,1,2)
    plot(stage_1_logAxis, stage_1_Log.stepSize)
    title('stage 1: step size')  
    % penalty parameter beta
    subplot(4,1,3)
    plot(stage_1_logAxis, stage_1_Log.beta)
    title('stage 1: penalty parameter beta')  
    % cost
    subplot(4,1,4)
    plot(stage_1_logAxis, log10(stage_1_Log.cost), 'g--',...
        stage_1_logAxis, log10(stage_1_Log.merit(:, 1)), 'k--',...
        stage_1_logAxis, log10(stage_1_Log.merit(:, 2)), 'r--')
    legend('cost', 'merit', 'merit(trial)')
    ylabel('log10')
    title('stage 1: cost and merit')        
    % KKT error
    figure(102)
    subplot(2,1,1)
    plot(stage_1_logAxis, log10(stage_1_Log.KKT_error(:, 1)), 'g--',...
        stage_1_logAxis, log10(stage_1_Log.KKT_error(:, 3)), 'k--',...
        stage_1_logAxis, log10(stage_1_Log.KKT_error(:, 5)), 'r--')
    legend('primal', 'dual(scaled)', 'complementary(scaled)')
    ylabel('log10')
    title('stage 1: KKT error')      
    subplot(2,1,2)
    plot(stage_1_logAxis, log10(stage_1_Log.FB_max(:, 1)), 'g--',...
        stage_1_logAxis, log10(stage_1_Log.FB_max(:, 2)), 'k--')
    legend('PSI_c', 'PSI_g')
    ylabel('log10')
    title('stage 1: FB max') 
end
% computation time
stage_1_Time = Info.Stage_1.Time;
figure(103)
bar(1, [stage_1_Time.gradEval, stage_1_Time.KKTEval, stage_1_Time.searchDirection, stage_1_Time.lineSearch, stage_1_Time.else, stage_1_Time.total])
legend({'gradEval', 'KKTEval', 'searchDirection', 'lineSearch', 'else', 'total'},'Location','northwest')
ylabel('time (s)'); title('stage 1: computation time')

%% stage 2: log and time
% log
Log = Info.Log;
logAxis = 0 : 1 : size(Log.cost, 1) - 1; % define log axis 
figure(105)
subplot(4,1,1)
plot(logAxis, log10(Log.param(:, 1)), 'g*-',...
    logAxis, log10(Log.param(:, 2)), 'bo-')
legend('s', '\sigma')
ylabel('log10')
title('parameter')
subplot(4,1,2)
plot(logAxis, Log.cost, 'r-')
title('cost')
subplot(4,1,3)
plot(logAxis, Log.VI_natural_residual, 'k-')
title('VI natural residual')
subplot(4,1,4)
stairs(logAxis(1 : end), [Log.timeElapsed(2 : end); Log.timeElapsed(end)], 'r-')
title('time Elapsed')

figure(106)
plot(logAxis, log10(Log.KKT_error(:, 1)), 'g--',...
    logAxis, log10(Log.KKT_error(:, 3)), 'k--',...
    logAxis, log10(Log.KKT_error(:, 5)), 'r--')
legend('primal', 'dual(scaled)', 'complementary(scaled)')
ylabel('log10')
title('KKT error')

% time
Time = Info.Time;
figure(103)
bar(1, [Time.non_interior_point, Time.Euler_predict_funcGradEval, Time.Euler_predict_stepEval,...
    Time.Newton_correct_funcGradEval, Time.Newton_correct_stepEval, Time.else, Time.total])
legend({'non interior point', 'Euler funcGradEval', 'Euler stepEval',...
    'Newton funcGradEval', 'Newton stepEval', 'else', 'total'},'Location','northwest')
ylabel('time (s)'); title('computation time')

end

