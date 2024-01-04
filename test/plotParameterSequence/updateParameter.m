clear all
clc

p_Init = [5e-1; 1e-1];
p_End = [1e-8; 1e-6];

kappa_times = 0.9;
kappa_exp = 1.1;

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

p = zeros(2, continuationStepNum + 1);
p_test = p_Init;
for i = 1 : continuationStepNum + 1
    p(:, i) = p_test;
    p_trial = min([kappa_times .* p_test, p_test.^kappa_exp], [], 2);
    p_test = max([p_trial, p_End], [], 2);    
end

figure(1)
semilogy(0:continuationStepNum , p(1, :), 'r*-',...
    0:continuationStepNum , p(2, :), 'go-')
grid on
legend('s', '\sigma')
xlabel('continuation step')
ylabel('parameter value')