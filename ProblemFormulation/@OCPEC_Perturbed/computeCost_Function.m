function L = computeCost_Function(OCPEC, Iterate, mode)
%computeCost_Function
%   Detailed explanation goes here

%%
nStages = OCPEC.nStages;
SC = OCPEC.StageCost;
TC = OCPEC.TerminalCost;
FRP = OCPEC.FRP;

%% compute cost function
L = zeros(1, nStages);
for n = 1 : nStages
    %
    tau_n = Iterate.tau(:, n);
    x_n = Iterate.x(:, n);
    p_n = Iterate.p(:, n);
    w_n = Iterate.w(:, n);
    switch mode
        case 'Regular'
            % compute stage cost 
            L_n = autoGen_LS(tau_n, x_n, p_n, w_n, SC.xRef(:, n), SC.tauRef(:, n), SC.xWeight, SC.tauWeight);
            % compute terminal cost
            if n == nStages
                L_T = autoGen_LT(tau_n, x_n, p_n, w_n, TC.xRef, TC.tauRef, TC.xWeight, TC.tauWeight);
                L_n = L_n + L_T;
            end
        case 'FRP'
            L_n = autoGen_LFRP(tau_n, x_n, p_n, w_n, FRP.ZRef(:, n), FRP.ZWeight(:, n));
    end
    % record
    L(:, n) = L_n;    
end

end

