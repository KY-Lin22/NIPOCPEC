function Option = createSolverOption()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% Option for stage 1: non-interior-point method
% initialization   
Option.Init.kappa_bound = 1e-2; % scaling parameter for single bound to evaluate perturbation (default 1e-2)
Option.Init.kappa_interval = 1e-2; % scaling parameter for interval to evaluate perturbation (default 1e-2)

% tolerance
Option.maxIterNum = 100;

Option.KKT_scaling_max = 1;
Option.tol.KKT_error_primal = 1e-6;
Option.tol.KKT_error_dual = 1e-4;
Option.tol.KKT_error_complementarity = 1e-4;
Option.tol.KKT_error_total = 1e-3;
Option.tol.dYNorm = 1e-6;

% singularity regularization parameter for KKT matrix 
Option.RegParam.nu_h = 1e-7; % rank-deficiency of equality constraint h
Option.RegParam.nu_c = 1e-8; % non-negative definite of diagonal matrix related to inequality constraint c
Option.RegParam.nu_g = 1e-8; % non-negative definite of diagonal matrix related to inequality constraint g
Option.RegParam.nu_H = 1e-8; % non-positive definite of Hessian matrix

% Lagrangian Hessian approximation
Option.KKT.Hessian_approximation = 'Gauss_Newton'; % 'Exact', 'Gauss_Newton'

% evaluate search direction

% merit line search
Option.LineSearch.betaInit = 1; % initial penalty parameter
Option.LineSearch.rho = 0.1; % desired extend for the negativity of merit function directional derivative
Option.LineSearch.stepSize_Min = 1e-4;
Option.LineSearch.stepSize_DecayRate = 0.5;% choose in (0,1)
Option.LineSearch.nu_D = 1e-4;% desired merit function reduction, default 1e-4

%% Option for stage 2: Euler-Newton continuation method
% init and end parameter
Option.Continuation.s_Init = 5e-1;
Option.Continuation.s_End = 1e-8;
Option.Continuation.sigma_Init = 1e-1;
Option.Continuation.sigma_End = 1e-6;
% update parameter 
Option.Continuation.kappa_times = 0.9;
Option.Continuation.kappa_exp = 1.1;
% tolerance
Option.Continuation.tol.KKT_error = 1e-4;
Option.Continuation.tol.VI_nat_res = 1e-2;
% Euler predict step

% Newton correction step
Option.Continuation.AdditionNewtonStep = 1;

end