function PHIvar = computeConstraint_Jacobian_PHI(OCPEC, Iterate, s)
%computeConstraint_Jacobian_PHI
%   Detailed explanation goes here
%%
nStages = OCPEC.nStages;
Dim = OCPEC.Dim;

%% compute constraint Jacobian
% define record
PHIvar = struct('PHItau', zeros(Dim.gamma, Dim.tau * nStages),...
                'PHIx', zeros(Dim.gamma, Dim.x * nStages),...
                'PHIp', zeros(Dim.gamma, Dim.p * nStages),...
                'PHIw', zeros(Dim.gamma, Dim.w * nStages));          
%  
for n = 1 : nStages
    p_n = Iterate.p(:, n);
    w_n = Iterate.w(:, n);
    
    [PHIp_n, PHIw_n] = autoGen_PHIp_PHIw(p_n, w_n, s);  
    
    PHIvar.PHIp(:, 1 + (n - 1) * Dim.p : n * Dim.p) = PHIp_n;
    PHIvar.PHIw(:, 1 + (n - 1) * Dim.w : n * Dim.w) = PHIw_n;     
end

end

