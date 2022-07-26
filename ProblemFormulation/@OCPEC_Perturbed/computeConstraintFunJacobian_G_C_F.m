function [Gvar, Cvar, Fvar] = computeConstraintFunJacobian_G_C_F(OCPEC, Iterate)
%computeConstraintFunJacobian_G_C_F
%   Detailed explanation goes here
%%
nStages = OCPEC.nStages;
Dim = OCPEC.Dim;

%% compute constraint Jacobian
% define record
Gvar = struct('Gtau', zeros(Dim.sigma, Dim.tau * nStages),...
              'Gx', zeros(Dim.sigma, Dim.x * nStages),...
              'Gp', zeros(Dim.sigma, Dim.p * nStages),...
              'Gw', zeros(Dim.sigma, Dim.w * nStages));
Cvar = struct('Ctau', zeros(Dim.eta, Dim.tau * nStages),...
              'Cx', zeros(Dim.eta, Dim.x * nStages),...
              'Cp', zeros(Dim.eta, Dim.p * nStages),...
              'Cw', zeros(Dim.eta, Dim.w * nStages));
Fvar = struct('Ftau', zeros(Dim.lambda, Dim.tau * nStages),...
              'Fx', zeros(Dim.lambda, Dim.x * nStages),...
              'Fp', zeros(Dim.lambda, Dim.p * nStages),...
              'Fw', zeros(Dim.lambda, Dim.w * nStages));     
%  
for n = 1 : nStages
    %
    tau_n = Iterate.tau(:, n);
    x_n = Iterate.x(:, n);
    p_n = Iterate.p(:, n);
    w_n = Iterate.w(:, n);
    
    [Gtau_n, Gx_n, Gp_n] = autoGen_Gtau_Gx_Gp(tau_n, x_n, p_n);    
    [Ctau_n, Cx_n, Cp_n, Cw_n] = autoGen_Ctau_Cx_Cp_Cw(tau_n, x_n, p_n, w_n);
    [ftau_n, fx_n, fp_n] = OCPEC.plant.computeStateEquationJacobian(tau_n, x_n, p_n); 
    % 
    Gvar.Gtau(:, 1 + (n - 1) * Dim.tau : n * Dim.tau) = Gtau_n;
    Gvar.Gx(:, 1 + (n - 1) * Dim.x : n * Dim.x) = Gx_n;
    Gvar.Gp(:, 1 + (n - 1) * Dim.p : n * Dim.p) = Gp_n;
    
    Cvar.Ctau(:, 1 + (n - 1) * Dim.tau : n * Dim.tau) = Ctau_n;
    Cvar.Cx(:, 1 + (n - 1) * Dim.x : n * Dim.x) = Cx_n;
    Cvar.Cp(:, 1 + (n - 1) * Dim.p : n * Dim.p) = Cp_n;
    Cvar.Cw(:, 1 + (n - 1) * Dim.w : n * Dim.w) = Cw_n;    
    
    Fvar.Ftau(:, 1 + (n - 1) * Dim.tau : n * Dim.tau) = ftau_n * OCPEC.timeStep;
    Fvar.Fx(:, 1 + (n - 1) * Dim.x : n * Dim.x) = fx_n * OCPEC.timeStep - eye(Dim.x);
    Fvar.Fp(:, 1 + (n - 1) * Dim.p : n * Dim.p) = fp_n * OCPEC.timeStep;    
      
end

end

