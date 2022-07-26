function Info = solutionExaminer(solver, solution, Record)
%solutionExaminer
%   Detailed explanation goes here
OCPEC = solver.OCPEC;
Option = solver.Option;

VI_mode = OCPEC.VI_mode;
plant = OCPEC.plant;
nStages = OCPEC.nStages;
Dim = OCPEC.Dim;

Tol = Option.Tolerance;
sEnd = Option.sEnd;
zEnd = Option.zEnd;
printLevel = Option.printLevel;

% Function evaluation
FunEval = solver.FunctionEvaluation(solution, sEnd, zEnd, 'Regular');
K = zeros(Dim.w, nStages);
for n = 1 : nStages
    K(:, n) = plant.computeVIFunc(solution.tau(:, n), solution.x(:, n), solution.p(:, n));
end

% Jacobian evaluation
FunEval.Lvar = OCPEC.computeCostFunJacobian(solution, 'Regular');
[FunEval.Gvar, FunEval.Cvar, FunEval.Fvar] = OCPEC.computeConstraintFunJacobian_G_C_F(solution);
FunEval.PHIvar = OCPEC.computeConstraintFunJacobian_PHI(solution, sEnd);

%% Solution Examiner: Cost, KKT Error
% total cost and KKT error
totalCost = sum(FunEval.L);
KKT_Residual = solver.computeKKT_Residual(solution, FunEval);
KKT_Error = solver.computeKKT_Error(solution, KKT_Residual);

%% Solution Examiner: Max Complementarity Residual between inequality constraints and its dual variables
% G and sigma
MaxCompResi_G_sigma = zeros(Dim.sigma, nStages); 
for n = 1 : nStages
    for i = 1 : Dim.sigma
        sigma_Vio = max(0, -solution.sigma(i, n));
        G_VioScale = min(1, max(0, solution.sigma(i, n)));
        MaxCompResi_G_sigma(i, n) = max(sigma_Vio, G_VioScale * max(FunEval.G(i, n), 0));
    end
end
MaxCompResi_G_sigma = norm(reshape(MaxCompResi_G_sigma, [], 1), Inf);
% PHI and gamma (only check in Reg_NCPs and Reg_Scholtes case)
MaxCompResi_PHI_gamma = zeros(Dim.gamma, nStages); 
if (strcmp(VI_mode,'Reg_NCPs')) || (strcmp(VI_mode, 'Reg_Scholtes'))
    for n = 1 : nStages
        for i = 1 : Dim.gamma
            gamma_Vio = max(0, -solution.gamma(i, n));
            PHI_VioScale = min(1, max(0, solution.gamma(i, n)));
            MaxCompResi_PHI_gamma(i, n) = max(gamma_Vio, PHI_VioScale * max(FunEval.PHI(i, n), 0));
        end
    end
end
MaxCompResi_PHI_gamma = norm(reshape(MaxCompResi_PHI_gamma, [], 1), Inf);

%% Solution Examiner: Max Constraint Residual/Violation
% Inequality-Type Constraint Violation
G_residual = min([zeros(Dim.sigma * nStages, 1), reshape(FunEval.G, [], 1)], [], 2);% check G >= 0
r_ineq_G = norm(G_residual, Inf);

% Equality-Type Constraint Residual
r_eq_C = norm(reshape(FunEval.C, [], 1), Inf); % check C = 0
r_eq_F = norm(reshape(FunEval.F, [], 1), Inf); % check F = 0

% Equilibrium Constraint Violation (inequality violation)
pK_residual = zeros(Dim.p, nStages); % for NCP check p >= 0 and K >=0; for VI check l <= p <= u
for n = 1 : nStages
    for i = 1 : Dim.p
        if (plant.l(i) == 0) && (plant.u(i) == Inf)
            % nonlinear complementary problem  
            p_Vio = abs(min(0, solution.p(i, n)));
            K_Vio = abs(min(0, K(i, n)));
            pK_residual(i, n) = max(p_Vio, K_Vio);                       
        else
            % box constraint variation inequality
            pl_Vio = abs(min(0, solution.p(i, n) - plant.l(i)));
            up_Vio = abs(min(0, plant.u(i) - solution.p(i, n)));                       
            pK_residual(i, n) = max(pl_Vio, up_Vio);
        end
    end
end
r_eqlb_ineq= norm(reshape(pK_residual, [], 1), Inf);

% Equilibrium Constraint Violation (complementarity violation)
pK_comp = zeros(Dim.p, nStages);
for n = 1 : nStages
    for i = 1 : Dim.p
        l_Vio = max(0, plant.l(i) - solution.p(i, n));
        K_l_VioScale = min(1, max(0, solution.p(i, n) - plant.l(i)));
        l_ComlVio = max(l_Vio, K_l_VioScale * max(K(i, n), 0));
        u_Vio = max(0, solution.p(i, n) - plant.u(i));
        K_u_VioScale = min(1, max(0, plant.u(i) - solution.p(i, n)));
        u_ComlVio = max(u_Vio, K_u_VioScale * max(-K(i, n), 0));
        pK_comp(i, n) = max(l_ComlVio, u_ComlVio);
    end
end
r_eqlb_comp = norm(reshape(pK_comp, [], 1), Inf);

%% creat Info and print message
% creat Info
Info.iterProcess = Record;
Info.solutionMsg.totalCost = totalCost;
Info.solutionMsg.KKT_Error = KKT_Error;
Info.solutionMsg.MaxCompResi_G_sigma = MaxCompResi_G_sigma;
Info.solutionMsg.MaxCompResi_PHI_gamma = MaxCompResi_PHI_gamma;
Info.solutionMsg.r_ineq_G = r_ineq_G;
Info.solutionMsg.r_eq_C = r_eq_C;
Info.solutionMsg.r_eq_F = r_eq_F;
Info.solutionMsg.r_eqlb_ineq = r_eqlb_ineq;
Info.solutionMsg.r_eqlb_comp = r_eqlb_comp;

% print message
disp('Done!');
disp('*------------------ Solution Information ------------------*')
% 1 terminal condition message
switch Record.terminalStatus
    case 1
        disp('1. Terminal Status: Solver finds the optimal solution')
        if (printLevel == 1) || (printLevel == 2) 
            if Record.terminalCond.sz
                disp('- This optimal solution is obtained with the desired perturbed parameters')
            end
            if Record.terminalCond.KKT_T
                disp('- This optimal solution satisfies the desired tolerance of total KKT error')
            end
            if Record.terminalCond.KKT_F
                disp('- This optimal solution satisfies the desired tolerance of feasibility')
            end
            if Record.terminalCond.KKT_S
                disp('- This optimal solution satisfies the desired tolerance of stationarity')
            end
        end
        
    case 0
        disp('1. Terminal Status: Solver fails to find the optimal solution')
        if (printLevel == 1) || (printLevel == 2)
            if Record.terminalCond.maxIterNum
                disp('- Failure Reasons: Solver has exceeded the maximum number of iterations as specified by the option')
            end
            if Record.terminalCond.LSwithoutFRP
                disp('- Failure Reasons: Line Search fails to find a new iterate with acceptable stepsize, while Feasibility Restoration Phase does not be employed')
            end
            if Record.terminalCond.LSwithFRP
                disp('- Failure Reasons: Feasibility Restoration Phase fails to find a less infeasibility solution')
            end
        end
end

% 2 Iteration Process Message
disp('2. Iteration Process Message')
disp(['- Iterations: ................', num2str(Record.iterNum)])
disp(['- TimeElapsed: ...............', num2str(Record.Time.total,'%10.3f'), 's'])
disp(['- AverageTime: ...............', num2str(1000 * Record.Time.total /Record.iterNum,'%10.2f'), ' ms/Iter'])

% 3 Solution Message
if (printLevel == 1) || (printLevel == 2)
    disp('3. Solution Message')
    disp('(1) Cost and KKT')
    disp(['- Total Cost: ................', num2str(totalCost,'%10.3f'), '; '])
    disp(['- KKT Error(Total): ..........', num2str(KKT_Error.Total,'%10.3e'),        ' / ', num2str(Tol.KKT_Error_Total,'%10.3e'), ' (opt / tol)'])
    disp(['           (Feasibility): ....', num2str(KKT_Error.Feasibility,'%10.3e'),  ' / ', num2str(Tol.KKT_Error_Feasibility,'%10.3e'), ' (opt / tol)'])
    disp(['           (Stationarity): ...', num2str(KKT_Error.Stationarity,'%10.3e'), ' / ', num2str(Tol.KKT_Error_Stationarity,'%10.3e'), ' (opt / tol)'])
    
    disp('(2) Max Complementarity Residual between inequality constraints and its dual variables')
    disp(['- G * sigma: .................',   num2str(MaxCompResi_G_sigma,'%10.3e')])
    disp(['- PHI * gamma: ...............', num2str(MaxCompResi_PHI_gamma,'%10.3e')])
    
    disp('(3) Constraint Residual/Violation')
    disp(['- Inequality-Type: ...........',   num2str(r_ineq_G,'%10.3e'), '(G >= 0); ']);
    disp(['- Equality-Type: .............',   num2str(r_eq_C,'%10.3e'), '(C = 0); ']);
    disp(['                 .............',   num2str(r_eq_F,'%10.3e'), '(F = 0); ']);
    disp(['- Equilibrium(ineq vio): .....',   num2str(r_eqlb_ineq,'%10.3e'), '(NCP: p >= 0, K >= 0; VI:l <= p <= u); '])
    disp(['             (comp vio): .....',   num2str(r_eqlb_comp,'%10.3e'), '(p * K = 0); '])
end
disp('*----------------------------------------------------------*')

end

