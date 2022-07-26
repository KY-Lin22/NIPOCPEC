function [dY, Info] = SearchDirection_Riccati(solver, KKT_Residual, KKT_Matrix)
%SearchDirection_Riccati
%   Detailed explanation goes here

TimeStart = tic;

OCPEC = solver.OCPEC;
Option = solver.Option;
nStages = OCPEC.nStages;
Dim = OCPEC.Dim;

% regroup KKT residual
T = [KKT_Residual.G_Fsb;...
     KKT_Residual.C_Fsb;...
     KKT_Residual.F_Fsb;...
     KKT_Residual.PHI_Fsb;...
     KKT_Residual.HtauT;...
     KKT_Residual.HxTlambdaNext;...
     KKT_Residual.HpT;...
     KKT_Residual.HwT];
%
J  = KKT_Matrix.J;
BL = KKT_Matrix.BL; 
BU = KKT_Matrix.BU;

%% backward recursion for Frak_A and Frak_b
Frak_A = zeros(Dim.Y, Dim.Y * nStages);
Frak_b = zeros(Dim.Y, nStages);
% initialization
J_N = J(:, Dim.Y * (nStages - 1) + 1: Dim.Y * nStages);
T_N = T(:, nStages);
Frak_A(:, Dim.Y * (nStages - 1) + 1: Dim.Y * nStages) = J_N;
Frak_b(:, nStages) = T_N;
for n = nStages - 1 : -1 : 1
    % load 
    J_n = J(:, Dim.Y * (n - 1) + 1 : Dim.Y * n);
    T_n = T(:, n);
    Frak_A_nNext = Frak_A(:, Dim.Y * n + 1 : Dim.Y * (n + 1));
    Frak_b_nNext = Frak_b(:, n + 1);
    % compute Frak_A_n and Frak_b_n
    if Option.employPinv
        Frak_A_n = J_n - BU * pinv(Frak_A_nNext) * BL;
        Frak_b_n = T_n - BU * pinv(Frak_A_nNext) * Frak_b_nNext;
    else
        Frak_A_n = J_n - BU * (Frak_A_nNext\BL);
        Frak_b_n = T_n - BU * (Frak_A_nNext\Frak_b_nNext);
    end
    % record
    Frak_A(:, Dim.Y * (n - 1) + 1 : Dim.Y * n) = Frak_A_n;
    Frak_b(:, n) = Frak_b_n;
end

%% forward recursion for dY
dY = zeros(Dim.Y, nStages);
dY_nPrev = zeros(Dim.Y, 1);
for n = 1 : nStages
    %
    Frak_A_n = Frak_A(:, Dim.Y * (n - 1) + 1 : Dim.Y * n);
    Frak_b_n = Frak_b(:, n);
    % compute search direction dY_n
    if Option.employPinv
        dY(:, n) = - pinv(Frak_A_n) * (Frak_b_n + BL * dY_nPrev);
    else
        dY(:, n) = - Frak_A_n\(Frak_b_n + BL * dY_nPrev);
    end 
    dY_nPrev = dY(:, n);
end

%% 
TimeElapsed = toc(TimeStart);
Info.Time = TimeElapsed;

end

