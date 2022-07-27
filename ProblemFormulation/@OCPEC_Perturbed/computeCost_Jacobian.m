function Lvar = computeCost_Jacobian(OCPEC, Iterate, mode)
%computeCost_Jacobian
%   Detailed explanation goes here

%%
nStages = OCPEC.nStages;
Dim = OCPEC.Dim;
SC = OCPEC.StageCost;
TC = OCPEC.TerminalCost;
FRP = OCPEC.FRP;

%% compute cost function Jacobian
% define record
Lvar = struct('Ltau', zeros(1, Dim.tau * nStages),...
              'Lx', zeros(1, Dim.x * nStages),...
              'Lp', zeros(1, Dim.p * nStages),...
              'Lw', zeros(1, Dim.w * nStages));
% compute cost Jacobian
for n = 1 : nStages
    %
    tau_n = Iterate.tau(:, n);
    x_n = Iterate.x(:, n);
    p_n = Iterate.p(:, n); 
    w_n = Iterate.w(:, n);
    switch mode
        case 'Regular'
            % compute stage cost Jacobian
            [Ltau_n, Lx_n, Lp_n, Lw_n] =...
                autoGen_LStau_LSx_LSp_LSw(tau_n, x_n, p_n, w_n, SC.xRef(:, n), SC.tauRef(:, n), SC.xWeight, SC.tauWeight);
            % compute terminal cost Jacobian
            if n == nStages
                [LTtau, LTx, LTp, LTw] =...
                    autoGen_LTtau_LTx_LTp_LTw(tau_n, x_n, p_n, w_n, TC.xRef, TC.tauRef, TC.xWeight, TC.tauWeight);
                Ltau_n = Ltau_n + LTtau;
                Lx_n = Lx_n + LTx;
                Lp_n = Lp_n + LTp;
                Lw_n = Lw_n + LTw;
            end
        case 'FRP'
            [Ltau_n, Lx_n, Lp_n, Lw_n] =...
                autoGen_LFRPtau_LFRPx_LFRPp_LFRPw(tau_n, x_n, p_n, w_n, FRP.ZRef(:, n),  FRP.ZWeight(:, n)); 
    end      
    % record
    Lvar.Ltau(:, 1 + (n - 1) * Dim.tau : n * Dim.tau) = Ltau_n;
    Lvar.Lx(:, 1 + (n - 1) * Dim.x : n * Dim.x) = Lx_n;
    Lvar.Lp(:, 1 + (n - 1) * Dim.p : n * Dim.p) = Lp_n;
    Lvar.Lw(:, 1 + (n - 1) * Dim.w : n * Dim.w) = Lw_n;  
end

end

