function value = computeFB_minusInvPSIbPSI(solver, PSIb, PSI)
%computeFB_minusInvPSIbPSI
%   Detailed explanation goes here

dim_PSI = size(PSI, 1);
nStages = size(PSI, 2);
nu_G = solver.Option.RegularParam.nu_G;

value = zeros(dim_PSI, nStages);
for n = 1 : nStages
    PSI_n = PSI(:, n);
    PSIb_n = PSIb(:, (n - 1) * dim_PSI + 1 : n * dim_PSI);
%     value(:, n) = - PSIb_n\PSI_n;   
    value(:, n) = - ( (diag(PSIb_n) - nu_G * ones(dim_PSI, 1) ).\ PSI_n);
end

end

