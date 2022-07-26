function showInfo(solver)
% showInfo
%   show information of dedicated NIP OCPEC solver
Option = solver.Option;

disp('*------------------------ NIPOCPEC Solver Information -----------------------*')
disp('1. Basic Options')
disp(['- maxIterNum: .............................', num2str(Option.maxIterNum)])
disp(['- KKT Error Tolerance(Total): .............', num2str(Option.Tolerance.KKT_Error_Total)])
disp(['                     (Feasibility): .......', num2str(Option.Tolerance.KKT_Error_Feasibility)])
disp(['                     (Stationarity): ......', num2str(Option.Tolerance.KKT_Error_Stationarity)])
disp('2. Options for Function and Jacobian Evaluation')
disp(['- Hessian Approximation Method: ...........', Option.HessianApproximation])
disp(['- Singularity Regular Parameter(J): .......', num2str(Option.RegularParam.nu_J)])
disp(['                               (G): .......', num2str(Option.RegularParam.nu_G)])
disp(['                               (H): .......', num2str(Option.RegularParam.nu_H)])
disp('3. Options for Search Direction Evaluation')
disp(['- Employ pinv .............................', mat2str(Option.employPinv)])
disp('4. Options for Line Search')
disp(['- Employ Second Order Correction: .........', mat2str(Option.employSecondOrderCorrection)]);
disp(['- stepSize_Min: ...........................', num2str(Option.LineSearch.stepSize_Min)])
disp('5. Options for Feasibility Restoration Phase')
disp(['- Employ Feasibility Restoration Phase: ...', mat2str(Option.employFeasibilityRestorationPhase)])
disp(['- maxIterNum: .............................', num2str(Option.FRP.maxIterNum)]);
disp(['- stepSize_Min: ...........................', num2str(Option.FRP.stepSize_Min)]);
disp('6. Options for Perturbed Parameter')
disp(['- sInit: ..................................', num2str(Option.sInit)])
disp(['- sEnd: ...................................', num2str(Option.sEnd)])
disp(['- zInit: ..................................', num2str(Option.zInit)])
disp(['- zEnd: ...................................', num2str(Option.zEnd)])
disp('*----------------------------------------------------------------------------*')
end

