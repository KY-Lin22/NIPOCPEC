function KKT_Matrix = computeKKT_Matrix(solver, FunEval)
%computeKKT_Matrix
%   Detailed explanation goes here

OCPEC = solver.OCPEC;
Option = solver.Option;
Dim = OCPEC.Dim;
nStages = OCPEC.nStages;
% regular parameter
nu_G = Option.RegularParam.nu_G;
nu_J = Option.RegularParam.nu_J;
nu_H = Option.RegularParam.nu_H;
% init
KKT_Matrix = struct('J', [], 'BL', [], 'BU', []);

%% KKT diagnal matrix: J
KKT_Matrix.J = zeros(Dim.Y, Dim.Y * nStages); 
for n = 1 : nStages
    % load Hessian and Constraint Jacobian
    Gvar_n = [FunEval.Gvar.Gtau(:, (n - 1) * Dim.tau + 1 : n * Dim.tau),...
        FunEval.Gvar.Gx(:, (n - 1) * Dim.x + 1 : n * Dim.x),...
        FunEval.Gvar.Gp(:, (n - 1) * Dim.p + 1 : n * Dim.p),...
        zeros(Dim.sigma, Dim.w)];  
    Cvar_n = [FunEval.Cvar.Ctau(:, (n - 1) * Dim.tau + 1 : n * Dim.tau),...
        FunEval.Cvar.Cx(:, (n - 1) * Dim.x + 1 : n * Dim.x),...
        FunEval.Cvar.Cp(:, (n - 1) * Dim.p + 1 : n * Dim.p),...
        FunEval.Cvar.Cw(:, (n - 1) * Dim.w + 1 : n * Dim.w)];    
    Fvar_n = [FunEval.Fvar.Ftau(:, (n - 1) * Dim.tau + 1 : n * Dim.tau),...
        FunEval.Fvar.Fx(:, (n - 1) * Dim.x + 1 : n * Dim.x),...
        FunEval.Fvar.Fp(:, (n - 1) * Dim.p + 1 : n * Dim.p),...
        zeros(Dim.lambda, Dim.w)];     
    PHIvar_n = [zeros(Dim.gamma, Dim.tau + Dim.x),...
        FunEval.PHIvar.PHIp(:, (n - 1) * Dim.p + 1 : n * Dim.p),...
        FunEval.PHIvar.PHIw(:, (n - 1) * Dim.w + 1 : n * Dim.w);];     
    Hessian_n = FunEval.Hessian(:, (n - 1) * Dim.Z + 1 : n * Dim.Z);     
    
    % D_n
    PSIgSigma_n = FunEval.PSIgSigma(:, (n - 1) * Dim.sigma + 1 : n * Dim.sigma);
    PSIgG_n     = FunEval.PSIgG(:, (n - 1) * Dim.sigma + 1 : n * Dim.sigma);
%     D_n = - PSIgG_n\PSIgSigma_n - nu_G * eye(Dim.sigma);  
    d_n = - ( ( diag(PSIgG_n) - nu_G * ones(Dim.sigma, 1) ).\ ( diag(PSIgSigma_n) - nu_G * ones(Dim.sigma, 1) ));
    D_n = diag(d_n);
    % J_n
    if (strcmp(OCPEC.VI_mode,'Reg_NCPs')) || (strcmp(OCPEC.VI_mode, 'Reg_Scholtes'))
        PSIphiGamma_n = FunEval.PSIphiGamma(:, (n - 1) * Dim.gamma + 1 : n * Dim.gamma);
        PSIphiPHI_n   = FunEval.PSIphiPHI(:, (n - 1) * Dim.gamma + 1 : n * Dim.gamma);
%         E_n = - PSIphiPHI_n\PSIphiGamma_n  - nu_G * eye(Dim.gamma);
        e_n = - ( (diag(PSIphiPHI_n) - nu_G * ones(Dim.gamma, 1)) .\ ( diag(PSIphiGamma_n) - nu_G * ones(Dim.gamma, 1) ));
        E_n = diag(e_n);
        dim_EtLam = Dim.eta + Dim.lambda;
        Avar_n = [Cvar_n; Fvar_n];
        J_n = [D_n,                          zeros(Dim.sigma, dim_EtLam),  zeros(Dim.sigma, Dim.gamma),  -Gvar_n;...
               zeros(dim_EtLam, Dim.sigma),  -nu_J * eye(dim_EtLam),       zeros(dim_EtLam, Dim.gamma),  Avar_n;...
               zeros(Dim.gamma, Dim.sigma),  zeros(Dim.gamma, dim_EtLam),  E_n,                          -PHIvar_n;...
               -Gvar_n',                     Avar_n',                      -PHIvar_n',                   Hessian_n + nu_H * eye(Dim.Z)];       
    elseif strcmp(OCPEC.VI_mode,  'SmoothingEquation')
        dim_EtLamGam = Dim.eta + Dim.lambda + Dim.gamma;
        Avar_n = [Cvar_n; Fvar_n; PHIvar_n];
        J_n = [D_n,                            zeros(Dim.sigma, dim_EtLamGam), -Gvar_n;...
               zeros(dim_EtLamGam, Dim.sigma), -nu_J * eye(dim_EtLamGam),      Avar_n;...
               -Gvar_n',                       Avar_n',                        Hessian_n + nu_H * eye(Dim.Z)];
    end
    
    KKT_Matrix.J(:, Dim.Y * (n - 1) + 1 : Dim.Y * n) = J_n;
end

%% KKT off-diagnal matrix: BL and BU
BL = [zeros(Dim.sigma + Dim.eta, Dim.Y);...
    zeros(Dim.lambda, Dim.LAMBDA + Dim.tau), eye(Dim.x), zeros(Dim.lambda, Dim.p + Dim.w);...
    zeros(Dim.gamma + Dim.Z, Dim.Y)];
KKT_Matrix.BL = BL;  
KKT_Matrix.BU = BL';

end
