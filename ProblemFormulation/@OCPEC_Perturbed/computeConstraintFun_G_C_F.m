function [G, C, F] = computeConstraintFun_G_C_F(OCPEC, Iterate)
%computeConstraintFun_G_C_F
%   Detailed explanation goes here

nStages = OCPEC.nStages;
Dim = OCPEC.Dim;
xPrev = [OCPEC.InitState, Iterate.x(:, 1 : end - 1)];

%% compute constraint function
G = zeros(Dim.sigma, nStages);
C = zeros(Dim.eta, nStages);
F = zeros(Dim.lambda, nStages);
for n = 1 : nStages
    %
    xPrev_n = xPrev(:, n);
    tau_n = Iterate.tau(:, n);
    x_n = Iterate.x(:, n);
    p_n = Iterate.p(:, n);
    w_n = Iterate.w(:, n);
    f_n = OCPEC.plant.computeStateEquation(tau_n, x_n, p_n);   
    %
    G(:, n) = autoGen_G(tau_n, x_n, p_n);
    C(:, n) = autoGen_C(tau_n, x_n, p_n, w_n);
    F(:, n) = xPrev_n - x_n + f_n * OCPEC.timeStep;
end

end

