function [P, l_Max] = createParameterSequence(self)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% load option parameter 
s_Init = self.Option.Continuation.s_Init;
s_End = self.Option.Continuation.s_End;
sigma_Init = self.Option.Continuation.sigma_Init;
sigma_End = self.Option.Continuation.sigma_End;
kappa_times = self.Option.Continuation.kappa_times;
kappa_exp = self.Option.Continuation.kappa_exp;

% formulate init and end parameter vector
p_Init = [s_Init; sigma_Init];
p_End = [s_End; sigma_End];

% evaluate parameter sequence and the max number of continuation step
p_l = p_Init;
P = p_l;
l_Max = 0;
while true
    if all(p_l == p_End)
        break
    else
        p_trial = min([kappa_times.*p_l, p_l.^kappa_exp], [], 2);
        p_l = max([p_trial, p_End], [], 2);
        P = [P, p_l];
        l_Max = l_Max + 1;
    end
end

end