function Sensitivity_Matrix = evaluateSensitivity_Matrix(self, PSI_c_grad_sigma,...
    PSI_g_grad_ineq, PSI_g_grad_sigma)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%   Sensitivity_Matrix has the form: 
% [zeros(Dim.h, 1),  zeros(Dim.h, 1);...
%  zeros(Dim.c, 1),  PSI_c_grad_sigma;...
%  PSI_g_grad_s,     PSI_g_grad_sigma;...
%  zeros(Dim.z, 1),  zeros(Dim.z, 1)]
%
Y_Node = cumsum([self.NLP.Dim.h, self.NLP.Dim.c, self.NLP.Dim.g, self.NLP.Dim.z]);
% PSI_g_grad_s
PSI_g_grad_s = PSI_g_grad_ineq * ones(self.NLP.Dim.g, 1);

%% extract nonzero index and value
% PSI_c_grad_sigma
[i_PSI_c_grad_sigma, j_PSI_c_grad_sigma, s_PSI_c_grad_sigma] = find(PSI_c_grad_sigma);
i_PSI_c_grad_sigma = i_PSI_c_grad_sigma + Y_Node(1);
j_PSI_c_grad_sigma = j_PSI_c_grad_sigma + 1;

% PSI_g_grad_s
[i_PSI_g_grad_s, j_PSI_g_grad_s, s_PSI_g_grad_s] = find(PSI_g_grad_s);
i_PSI_g_grad_s = i_PSI_g_grad_s + Y_Node(2);

% PSI_g_grad_sigma
[i_PSI_g_grad_sigma, j_PSI_g_grad_sigma, s_PSI_g_grad_sigma] = find(PSI_g_grad_sigma);
i_PSI_g_grad_sigma = i_PSI_g_grad_sigma + Y_Node(2);
j_PSI_g_grad_sigma = j_PSI_g_grad_sigma + 1;

%% assemble Sensitivity matrix
i_Sensitivity = [i_PSI_c_grad_sigma; i_PSI_g_grad_s; i_PSI_g_grad_sigma];
j_Sensitivity = [j_PSI_c_grad_sigma; j_PSI_g_grad_s; j_PSI_g_grad_sigma];
s_Sensitivity = [s_PSI_c_grad_sigma; s_PSI_g_grad_s; s_PSI_g_grad_sigma];

Sensitivity_Matrix = sparse(i_Sensitivity, j_Sensitivity, s_Sensitivity,...
    Y_Node(4), 2, length(s_Sensitivity));
end

