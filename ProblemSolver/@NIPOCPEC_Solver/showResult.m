function showResult(solver, Info)
%showResult
%   Detailed explanation goes here

Option= solver.Option;
iterProcess = Info.iterProcess;
iterNum = iterProcess.iterNum;
iterCounterAxis = 0 : 1 : iterProcess.iterNum;% define iterCounter axis

%% plot iteration process
figure(31)
% total cost
subplot(5,1,1)
plot(iterCounterAxis, iterProcess.totalCost(1 : iterNum + 1, 1))
title('total cost')

% KKT error
subplot(5,1,2)
plot(iterCounterAxis, log10(iterProcess.KKT_Error.Total(1 : iterNum + 1, 1)), 'k-')
hold on
plot(iterCounterAxis, log10(iterProcess.KKT_Error.Feasibility(1 : iterNum + 1, 1)), 'r--')
hold on
plot(iterCounterAxis, log10(iterProcess.KKT_Error.Stationarity(1 : iterNum + 1, 1)), 'g:')
legend('Total', 'Feasibility', 'Stationarity')
ylabel('log10'); title('KKT Error')

% perturbed parameter
subplot(5,1,3)
plot(iterCounterAxis, iterProcess.s(1 : iterNum + 1, 1), 'k-')
hold on
plot(iterCounterAxis, iterProcess.z(1 : iterNum + 1, 1), 'b-')
legend('s', 'z')
title('perturbed parameter ')

% penalty parameter and merit
subplot(5,1,4)
plot(iterCounterAxis, log10(iterProcess.beta(1: iterNum + 1, 1)), 'k-')
hold on
plot(iterCounterAxis, log10(iterProcess.Merit(1 : iterNum + 1, 1)), 'r--')
legend('beta', 'merit')
ylabel('log10'); title('penalty parameter and merit')

% stepsize
subplot(5,1,5)

MLS_index = find(strcmp(iterProcess.IterateType, 'MLS'));
SOC_index = find(strcmp(iterProcess.IterateType, 'SOC'));
FRP_index = find(strcmp(iterProcess.IterateType, 'FRP'));

plot(iterCounterAxis, iterProcess.stepSize(1: iterNum + 1, 1), 'k-', 'Marker', '.', 'MarkerIndices', MLS_index)
hold on
plot(iterCounterAxis, iterProcess.stepSize(1: iterNum + 1, 1), 'k-', 'Marker', '*', 'MarkerIndices', SOC_index)
hold on
plot(iterCounterAxis, iterProcess.stepSize(1: iterNum + 1, 1), 'k-', 'Marker', 'o', 'MarkerIndices', FRP_index)
legend('MLS', 'SOC', 'FRP')
xlabel('Iteration Number'); 
title('step size')

% computation time 
figure(41)
Time = [iterProcess.Time.JacobianHessian;...
    iterProcess.Time.KKT;...
    iterProcess.Time.SearchDirection;...
    iterProcess.Time.LineSearch;...
    iterProcess.Time.FRP;...
    iterProcess.Time.else;...
    iterProcess.Time.total];
bar(1,Time)
legend({'JacobianHessian', 'KKT', 'Search direction', 'LineSearch', 'FRP', 'else', 'total'},'Location','northwest')
ylabel('computation time (s)'); title('computation time')

end

