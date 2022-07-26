function codeGen(solver)
%codeGen
%   generate code files about Hessian

disp('Generate Code Files about Hessian..');
Option = solver.Option;
OCPEC = solver.OCPEC;

%% define Hamiltion
switch Option.HessianApproximation
    case 'GaussNewton'
        H     = OCPEC.LS;
        H_FRP = OCPEC.LFRP;        
    case 'Newton'
        if plant.computeInvM
            % force GaussNewton
            H     = OCPEC.LS;
            H_FRP = OCPEC.LFRP;
        else
            if (strcmp(OCPEC.VI_mode, 'Reg_NCPs')) || (strcmp(OCPEC.VI_mode, 'Reg_Scholtes'))
                H     = OCPEC.LS    - OCPEC.sigma' * OCPEC.G + OCPEC.eta' * OCPEC.C...
                    + OCPEC.lambda' * (OCPEC.plant.f * OCPEC.timeStep - OCPEC.x) - OCPEC.gamma' * OCPEC.PHI;
                H_FRP = OCPEC.LFRP - OCPEC.sigma' * OCPEC.G + OCPEC.eta' * OCPEC.C...
                    + OCPEC.lambda' * (OCPEC.plant.f * OCPEC.timeStep - OCPEC.x) - OCPEC.gamma' * OCPEC.PHI;
            elseif (strcmp(OCPEC.VI_mode, 'SmoothingEquation'))
                H     = OCPEC.LS    - OCPEC.sigma' * OCPEC.G + OCPEC.eta' * OCPEC.C...
                    + OCPEC.lambda' * (OCPEC.plant.f * OCPEC.timeStep - OCPEC.x) + OCPEC.gamma' * OCPEC.PHI;
                H_FRP = OCPEC.LFRP - OCPEC.sigma' * OCPEC.G + OCPEC.eta' * OCPEC.C...
                    + OCPEC.lambda' * (OCPEC.plant.f * OCPEC.timeStep - OCPEC.x) + OCPEC.gamma' * OCPEC.PHI;
            end
        end       
    otherwise
        % force 'GaussNewton'
        H     = OCPEC.LS;
        H_FRP = OCPEC.LFRP;        
end

%% compute Hessian of Hamiltion and LT
Hessian     = computeHessianSymFun(H, OCPEC);
Hessian_FRP = computeHessianSymFun(H_FRP, OCPEC);
Hessian_LT  = computeHessianSymFun(OCPEC.LT, OCPEC);

%% generate code files
% Hessian of Hamiltion
varsList = {OCPEC.sigma, OCPEC.eta, OCPEC.lambda, OCPEC.gamma,...
    OCPEC.tau, OCPEC.x, OCPEC.p, OCPEC.w,...
    OCPEC.symVar.xRef, OCPEC.symVar.tauRef, OCPEC.symVar.xWeight, OCPEC.symVar.tauWeight, OCPEC.symVar.s};
matlabFunction(Hessian,...
    'file','./autoGen_CodeFiles/autoGen_Hessian.m',...
    'vars', varsList,...
    'outputs',{'Hessian'},...
    'Optimize',OCPEC.codeOptimize);
% Hessian of LT
varsList_T = {OCPEC.sigma, OCPEC.eta, OCPEC.lambda, OCPEC.gamma,...
    OCPEC.tau, OCPEC.x, OCPEC.p, OCPEC.w,...
    OCPEC.symVar.xRef, OCPEC.symVar.tauRef, OCPEC.symVar.xWeight, OCPEC.symVar.tauWeight};
matlabFunction(Hessian_LT,...
    'file','./autoGen_CodeFiles/autoGen_Hessian_LT.m',...
    'vars', varsList_T,...
    'outputs',{'Hessian_LT'},...
    'Optimize',OCPEC.codeOptimize);
% Hessian of Hamiltion in FRP
varsList_FRP = {OCPEC.sigma, OCPEC.eta, OCPEC.lambda, OCPEC.gamma,...
    OCPEC.tau, OCPEC.x, OCPEC.p, OCPEC.w,...
    OCPEC.symVar.ZRef, OCPEC.symVar.ZWeight, OCPEC.symVar.s};
matlabFunction(Hessian_FRP,...
    'file','./autoGen_CodeFiles/autoGen_Hessian_FRP.m',...
    'vars', varsList_FRP,...
    'outputs',{'Hessian_FRP'},...
    'Optimize',OCPEC.codeOptimize);

disp('Done!')

end
%%
function Hessian = computeHessianSymFun(H, OCPEC)
% compute Jacobian
Htau = jacobian(H, OCPEC.tau);
Hx   = jacobian(H, OCPEC.x);
Hp   = jacobian(H, OCPEC.p);
Hw   = jacobian(H, OCPEC.w);
% compute Hessian
Htautau = jacobian(Htau, OCPEC.tau);
Htaux   = jacobian(Htau, OCPEC.x);
Htaup   = jacobian(Htau, OCPEC.p);
Htauw   = jacobian(Htau, OCPEC.w);

Hxx = jacobian(Hx, OCPEC.x);
Hxp = jacobian(Hx, OCPEC.p);
Hxw = jacobian(Hx, OCPEC.w);

Hpp = jacobian(Hp, OCPEC.p);
Hpw = jacobian(Hp, OCPEC.w);
    
Hww = jacobian(Hw, OCPEC.w);   

Hessian = [Htautau, Htaux, Htaup, Htauw;...
           Htaux',  Hxx,   Hxp,   Hxw;...
           Htaup',  Hxp',  Hpp,   Hpw;...
           Htauw',  Hxw',  Hpw',  Hww];
end
