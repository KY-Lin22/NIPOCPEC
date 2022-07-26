function codeGen(plant)
%codeGen
%   generate codes to compute dynamics

%% check symbolic function of dynamics
% state equation
if plant.computeInvM
    q_Dim = round(1/2 * plant.Dim.x);
    if ~all(size(plant.M) == [q_Dim, q_Dim])
    error('please specify plant.M with correct dimension')
    end
    if ~all(size(plant.H) == [q_Dim, 1])
    error('please specify plant.H with correct dimension')
    end    
else
    if ~all(size(plant.f) == [plant.Dim.x, 1])
    error('please specify plant.f with correct dimension')
    end    
end
% VI
if ~all(size(plant.K) == [plant.Dim.p, 1])
    error('please specify plant.K with correct dimension')
end
if ~all(size(plant.l) == [plant.Dim.p, 1])
    error('please specify plant.l with correct dimension')
end
if ~all(size(plant.u) == [plant.Dim.p, 1])
    error('please specify plant.u with correct dimension')
end

disp('Generate code files about dynamics...');

%% Generate code files about state equation
if plant.computeInvM
    % compute Jacobian 
    Mx = sym('Mx', [q_Dim, q_Dim * plant.Dim.x]);
    for i = 1 : plant.Dim.x
        Mx(:, (i - 1)*q_Dim + 1 : i * q_Dim) = formula(diff(plant.M, plant.x(i)));
    end
    Mx = symfun(Mx, plant.x);
    Htau = jacobian(plant.H, plant.tau);
    Hx = jacobian(plant.H, plant.x);
    Hp = jacobian(plant.H, plant.p);  
    
    matlabFunction(plant.M,...
        'file','./autoGen_CodeFiles/autoGen_M.m',...
        'vars',{plant.x},...
        'outputs',{'M'},...
        'Optimize',plant.codeOptimize);
    matlabFunction(plant.H,...
        'file','./autoGen_CodeFiles/autoGen_H.m',...
        'vars',{plant.tau, plant.x, plant.p},...
        'outputs',{'H'},...
        'Optimize',plant.codeOptimize);
    
    matlabFunction(Mx,...
        'file','./autoGen_CodeFiles/autoGen_Mx.m',...
        'vars',{plant.x},...
        'outputs',{'Mx'},...
        'Optimize',plant.codeOptimize);    
    matlabFunction(Htau, Hx, Hp,...
        'file','./autoGen_CodeFiles/autoGen_Htau_Hx_Hp.m',...
        'vars',{plant.tau, plant.x, plant.p},...
        'outputs',{'Htau', 'Hx', 'Hp'},...
        'Optimize',plant.codeOptimize);
else
    % compute jacobian
    ftau = jacobian(plant.f, plant.tau);
    fx   = jacobian(plant.f, plant.x);
    fp   = jacobian(plant.f, plant.p);

    matlabFunction(plant.f,...
        'file','./autoGen_CodeFiles/autoGen_f.m',...
        'vars',{plant.tau, plant.x, plant.p},...
        'outputs',{'f'},...
        'Optimize',plant.codeOptimize);    
    matlabFunction(ftau, fx, fp,...
        'file','./autoGen_CodeFiles/autoGen_ftau_fx_fp.m',...
        'vars',{plant.tau, plant.x, plant.p},...
        'outputs',{'ftau', 'fx', 'fp'},...
        'Optimize',plant.codeOptimize);    
end

%% Generate code files about VI
matlabFunction(plant.K,...
    'file','./autoGen_CodeFiles/autoGen_K.m',...
    'vars',{plant.tau, plant.x, plant.p},...
    'outputs',{'K'},...
    'Optimize',plant.codeOptimize);

disp('Done!')

end

