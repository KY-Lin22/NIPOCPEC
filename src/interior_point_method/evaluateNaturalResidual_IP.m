function natRes = evaluateNaturalResidual_IP(OCPEC, NLP, z_Opt)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*

%% extract solution
Z_Opt = reshape(z_Opt, NLP.Dim.z_Node(end), OCPEC.nStages);

X_Opt = Z_Opt(1 : NLP.Dim.z_Node(1), :);
U_Opt = Z_Opt(NLP.Dim.z_Node(1) + 1 : NLP.Dim.z_Node(2), :);
LAMBDA_Opt = Z_Opt(NLP.Dim.z_Node(2) + 1 : NLP.Dim.z_Node(3), :);

%% problem data for Euclidean projector
F_FuncObj_map = OCPEC.FuncObj.F.map(OCPEC.nStages);
F_Opt = full(F_FuncObj_map(X_Opt, U_Opt, LAMBDA_Opt));
w_Opt = LAMBDA_Opt - F_Opt;

%% evaluate Euclidean projector
switch OCPEC.VISetType
    case 'box_constraint'
        % construct a box-constraint projector to evaluate Euclidean projector
        w = SX.sym('w', OCPEC.Dim.lambda, 1); % projector input
        proj_w_func = max(OCPEC.bl, min(w, OCPEC.bu)); % projector function
        solver_singleStage = Function('EuclideanProjector', {w}, {proj_w_func}, {'w'}, {'proj_w_func'});
        EuclideanProjector = solver_singleStage.map(OCPEC.nStages);
        proj_w_Opt = full(EuclideanProjector(w_Opt));
    case 'nonnegative_orthant'
        % construct a max projector to evaluate Euclidean projector
        w = SX.sym('w', OCPEC.Dim.lambda, 1); % projector input
        proj_w_func = max(zeros(OCPEC.Dim.lambda, 1), w); % projector function
        solver_singleStage = Function('EuclideanProjector', {w}, {proj_w_func}, {'w'}, {'proj_w_func'});
        EuclideanProjector = solver_singleStage.map(OCPEC.nStages);
        proj_w_Opt = full(EuclideanProjector(w_Opt));        
end

%% evaluate natural residual
natRes = reshape(LAMBDA_Opt - proj_w_Opt, [], 1);
end


