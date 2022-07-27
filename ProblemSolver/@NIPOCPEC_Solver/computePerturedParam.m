function [s_k, z_k, FunEval_k] = computePerturedParam(solver, Iterate_k, FunEval_k, s, z)
%computePerturedParam
%   Detailed explanation goes here
OCPEC = solver.OCPEC;
Option = solver.Option;
sEnd = Option.sEnd;
zEnd = Option.zEnd;

kappa_s_times = 0.2;
kappa_s_exp = 1.5;
kappa_z_times = 0.2;
kappa_z_exp = 1.5;

% set scaling parameter and threshold 
kappa_F = 10; 
threshold_sz = max([s, z]);

% constraint violation (L Inf norm)
if (strcmp(OCPEC.VI_mode,'Reg_NCPs')) || (strcmp(OCPEC.VI_mode, 'Reg_Scholtes'))
    totalCstrVio_LInfNorm = norm(reshape([FunEval_k.PSIg; FunEval_k.C; FunEval_k.F; FunEval_k.PSIphi], [], 1), Inf); 
elseif strcmp(OCPEC.VI_mode, 'SmoothingEquation')
    totalCstrVio_LInfNorm = norm(reshape([FunEval_k.PSIg; FunEval_k.C; FunEval_k.F; FunEval_k.PHI], [], 1), Inf); 
end

%% update s and z
if totalCstrVio_LInfNorm < kappa_F * threshold_sz
    % close to the optimal solution, need to set them a smaller value
    s_k_trail = min([kappa_s_times .* s, s.^kappa_s_exp]);
    s_k = max([s_k_trail, sEnd]);
    z_k_trail = min([kappa_z_times .* z, z.^kappa_z_exp]);
    z_k = max([z_k_trail, zEnd]);
else
    % far way from the optimal solution, set them the previous value
    s_k = s;
    z_k = z;
end

%% update FunEval_k
if (s_k == s) && (z_k == z)
    % both s and z do not update, hence no function need to update
     
elseif (s_k ~= s) && (z_k == z)
    % only s is update
    FunEval_k.PHI = OCPEC.computeConstraint_Function_PHI(Iterate_k, s_k);
    if (strcmp(OCPEC.VI_mode,'Reg_NCPs')) || (strcmp(OCPEC.VI_mode, 'Reg_Scholtes'))
        [FunEval_k.PSIphi, FunEval_k.PSIphiGamma, FunEval_k.PSIphiPHI] =...
            solver.computeFB_Function_Jacobian(Iterate_k.gamma, FunEval_k.PHI, z_k);        
    end
    
elseif (s_k == s) && (z_k ~= z)
    % only z is update
    [FunEval_k.PSIg, FunEval_k.PSIgSigma, FunEval_k.PSIgG] =...
        solver.computeFB_Function_Jacobian(Iterate_k.sigma, FunEval_k.G, z_k);
    
    if (strcmp(OCPEC.VI_mode,'Reg_NCPs')) || (strcmp(OCPEC.VI_mode, 'Reg_Scholtes'))
        [FunEval_k.PSIphi, FunEval_k.PSIphiGamma, FunEval_k.PSIphiPHI] =...
            solver.computeFB_Function_Jacobian(Iterate_k.gamma, FunEval_k.PHI, z_k);
    end
    
else
    % both s and z update
    [FunEval_k.PSIg, FunEval_k.PSIgSigma, FunEval_k.PSIgG] =...
        solver.computeFB_Function_Jacobian(Iterate_k.sigma, FunEval_k.G, z_k);
    
    FunEval_k.PHI = OCPEC.computeConstraint_Function_PHI(Iterate_k, s_k);
    if (strcmp(OCPEC.VI_mode,'Reg_NCPs')) || (strcmp(OCPEC.VI_mode, 'Reg_Scholtes'))
        [FunEval_k.PSIphi, FunEval_k.PSIphiGamma, FunEval_k.PSIphiPHI] =...
            solver.computeFB_Function_Jacobian(Iterate_k.gamma, FunEval_k.PHI, z_k);
    end
        
end

end

