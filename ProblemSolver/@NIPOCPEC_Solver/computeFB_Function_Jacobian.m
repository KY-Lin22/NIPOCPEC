function [PSI, PSIa, PSIb] = computeFB_Function_Jacobian(solver, a, b, z)
%computeFB_Function_Jacobian
%   compute Fischer Burmeister Function and Jacobian

% check input
if ~all(size(a) == size(b))
    error('vector a and b should have the same dimension');
end

% function and jacobian computation
dim_PSI = size(a, 1);
nStages = size(a, 2);

PSI = zeros(dim_PSI, nStages);
PSIa = zeros(dim_PSI, dim_PSI * nStages);
PSIb = zeros(dim_PSI, dim_PSI * nStages);
for n = 1 : nStages
    a_n = a(:, n);
    b_n = b(:, n);
    PSI_n = zeros(dim_PSI, 1);
    PSIa_n = zeros(dim_PSI, dim_PSI);
    PSIb_n = zeros(dim_PSI, dim_PSI);
    for i = 1 : dim_PSI
        Sq = sqrt((a_n(i)).^2 + (b_n(i)).^2 + (z).^2);
        PSI_n(i) = Sq - a_n(i) - b_n(i);
        PSIa_n(i, i) = a_n(i)/Sq - 1;
        PSIb_n(i, i) = b_n(i)/Sq - 1;       
    end
    PSI(:, n) = PSI_n;
    PSIa(:, (n - 1) * dim_PSI + 1 : n * dim_PSI) = PSIa_n;
    PSIb(:, (n - 1) * dim_PSI + 1 : n * dim_PSI) = PSIb_n;
end

end
