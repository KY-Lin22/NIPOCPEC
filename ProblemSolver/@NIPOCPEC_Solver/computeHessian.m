function Hessian = computeHessian(solver, Iterate, s, mode)
%computeHessian
%   Detailed explanation goes here

OCPEC = solver.OCPEC;
nStages = OCPEC.nStages;
Dim = OCPEC.Dim;
SC = OCPEC.StageCost;
TC = OCPEC.TerminalCost;
FRP = OCPEC.FRP;

%% compute Hessian
Hessian = zeros(Dim.Z, Dim.Z * nStages);
for n = 1 : nStages
    % load
    sigma_n = Iterate.sigma(:, n);
    eta_n = Iterate.eta(:, n);
    lambda_n = Iterate.lambda(:, n);
    gamma_n = Iterate.gamma(:, n);   
    tau_n = Iterate.tau(:, n);
    x_n = Iterate.x(:, n);
    p_n = Iterate.p(:, n);
    w_n = Iterate.w(:, n);
    % compute
    switch mode
        case 'Regular'
            % compute Hessian of Hamiltion
            Hessian_n = autoGen_Hessian(sigma_n, eta_n, lambda_n, gamma_n, tau_n, x_n, p_n, w_n,...
                SC.xRef(:, n), SC.tauRef(:, n), SC.xWeight, SC.tauWeight, s);
            % compute Hessian of LT
            if n == nStages
                Hessian_LT = autoGen_Hessian_LT(sigma_n, eta_n, lambda_n, gamma_n, tau_n, x_n, p_n, w_n,...
                    TC.xRef, TC.tauRef, TC.xWeight, TC.tauWeight);
                Hessian_n = Hessian_n + Hessian_LT;
            end
        case 'FRP'
            Hessian_n = autoGen_Hessian_FRP(sigma_n, eta_n, lambda_n, gamma_n, tau_n, x_n, p_n, w_n,...
                FRP.ZRef(:, n), FRP.ZWeight(:, n), s);
    end
    % record
    Hessian(:, (n - 1) * Dim.Z + 1 : n * Dim.Z) = Hessian_n;        
end

end

