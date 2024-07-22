function [Y_l, Info] = EulerNewtonContinuationMethod(self, Y, p_l, p)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
timeStart = tic;

% Euler predict step
dp = p_l - p;
dY_Euler = full(self.FuncObj.dY_Euler(Y, p, dp));
Y_Euler = Y + dY_Euler;

% Newton corrector step
Y_Newton = Y_Euler;
for i = 1 : 1 + self.Option.Continuation.AdditionNewtonStep
    dY_Newton = full(self.FuncObj.dY_Newton(Y_Newton, p_l));
    Y_Newton = Y_Newton + dY_Newton;
end

% Info
timeElapsed = toc(timeStart);
terminal_status = 1;
terminal_msg = ['- Solver succeeds: ', 'because a new iterate found by Euler Newton continuation method'];
Y_l = Y_Newton;
Info.time = timeElapsed;
Info.terminal_status = terminal_status;
Info.terminal_msg = terminal_msg;

end