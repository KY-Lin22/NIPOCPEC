function PHI = computeConstraint_Function_PHI(OCPEC, Iterate, s)
%computeConstraint_Function_PHI
%   Detailed explanation goes here
nStages = OCPEC.nStages;
Dim = OCPEC.Dim;

%% compute constraint function
PHI = zeros(Dim.gamma, nStages);
for n = 1 : nStages
    p_n = Iterate.p(:, n);
    w_n = Iterate.w(:, n); 
    PHI(:, n) = autoGen_PHI(p_n, w_n, s);
end

end

