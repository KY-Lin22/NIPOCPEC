function KKT_Matrix = evaluateKKT_Matrix(self, h_grad, c_grad, g_grad, LAG_hessian,...
            PSI_c_grad_dual, PSI_c_grad_ineq, PSI_g_grad_dual, PSI_g_grad_ineq)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
%
%   KKT_Matrix has the form: 
%   [-nu_h * eye(Dim.h),    zeros(Dim.h, Dim.c),                   zeros(Dim.h, Dim.g),                   h_grad;... 
%    zeros(Dim.c, Dim.h),   PSI_c_grad_dual - nu_c * eye(Dim.c),   zeros(Dim.c, Dim.g),                   PSI_c_grad_z;... 
%    zeros(Dim.g, Dim.h),   zeros(Dim.g, Dim.c),                   PSI_g_grad_dual - nu_g * eye(Dim.g),   PSI_g_grad_z;...
%    h_grad',               -c_grad',                              -g_grad',                              Hessian + nu_H * eye(Dim.z)]
%
%%
% regularization parameter and Y node point
nu_h = self.Option.RegParam.nu_h;
nu_c = self.Option.RegParam.nu_c;
nu_g = self.Option.RegParam.nu_g;
nu_H = self.Option.RegParam.nu_H;
Y_Node = cumsum([self.NLP.Dim.h, self.NLP.Dim.c, self.NLP.Dim.g, self.NLP.Dim.z]);
% PSI_c_grad_z and PSI_g_grad_z
PSI_c_grad_z = PSI_c_grad_ineq * c_grad;
PSI_g_grad_z = PSI_g_grad_ineq * g_grad;
% diag vector of top left 3 X 3 block matrix
diagVec = [- nu_h * ones(self.NLP.Dim.h, 1);...
    diag(PSI_c_grad_dual) - nu_c * ones(self.NLP.Dim.c, 1);...
    diag(PSI_g_grad_dual) - nu_g * ones(self.NLP.Dim.g, 1)];

%% extract nonzero index and value
% top left 3 X 3 block matrix
i_diagVec = 1 : Y_Node(3);
j_diagVec = 1 : Y_Node(3);
s_diagVec = diagVec;

% h_grad
[i_h_grad, j_h_grad, s_h_grad] = find(h_grad);
j_h_grad = j_h_grad + Y_Node(3);

% PSI_c_grad_z
[i_PSI_c_grad_z, j_PSI_c_grad_z, s_PSI_c_grad_z] = find(PSI_c_grad_z);
i_PSI_c_grad_z = i_PSI_c_grad_z + Y_Node(1);
j_PSI_c_grad_z = j_PSI_c_grad_z + Y_Node(3);

% PSI_g_grad_z
[i_PSI_g_grad_z, j_PSI_g_grad_z, s_PSI_g_grad_z] = find(PSI_g_grad_z);
i_PSI_g_grad_z = i_PSI_g_grad_z + Y_Node(2);
j_PSI_g_grad_z = j_PSI_g_grad_z + Y_Node(3);

% h_grad'
[i_h_grad_T, j_h_grad_T, s_h_grad_T] = find(h_grad');
i_h_grad_T = i_h_grad_T + Y_Node(3);

% -c_grad'
[i_nega_c_grad_T, j_nega_c_grad_T, s_nega_c_grad_T] = find(-c_grad');
i_nega_c_grad_T = i_nega_c_grad_T + Y_Node(3);
j_nega_c_grad_T = j_nega_c_grad_T + Y_Node(1);

% -g_grad'
[i_nega_g_grad_T, j_nega_g_grad_T, s_nega_g_grad_T] = find(-g_grad');
i_nega_g_grad_T = i_nega_g_grad_T + Y_Node(3);
j_nega_g_grad_T = j_nega_g_grad_T + Y_Node(2);

% Hessian + nu_H * eye(Dim.z)
[i_hessian, j_hessian, s_hessian] = find(LAG_hessian + nu_H * speye(self.NLP.Dim.z));
i_hessian = i_hessian + Y_Node(3);
j_hessian = j_hessian + Y_Node(3);

%% assemble KKT matrix
i_KKT = [i_diagVec'; i_h_grad; i_PSI_c_grad_z; i_PSI_g_grad_z; i_h_grad_T; i_nega_c_grad_T; i_nega_g_grad_T; i_hessian];
j_KKT = [j_diagVec'; j_h_grad; j_PSI_c_grad_z; j_PSI_g_grad_z; j_h_grad_T; j_nega_c_grad_T; j_nega_g_grad_T; j_hessian];
s_KKT = [s_diagVec;  s_h_grad; s_PSI_c_grad_z; s_PSI_g_grad_z; s_h_grad_T; s_nega_c_grad_T; s_nega_g_grad_T; s_hessian];

KKT_Matrix = sparse(i_KKT, j_KKT, s_KKT, Y_Node(4), Y_Node(4), length(s_KKT));
end

