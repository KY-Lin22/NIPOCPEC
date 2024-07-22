function [KKT_error_primal, KKT_error_dual, KKT_error_complementary, KKT_error_total] = evaluateKKT_Error(self, Y, p)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here

% extract elements from Y and p
Y_Node = cumsum([self.NLP.Dim.h, self.NLP.Dim.c, self.NLP.Dim.g, self.NLP.Dim.z]);
gamma_h = Y(            1 : Y_Node(1), 1);
gamma_c = Y(Y_Node(1) + 1 : Y_Node(2), 1);
gamma_g = Y(Y_Node(2) + 1 : Y_Node(3), 1);
z       = Y(Y_Node(3) + 1 : Y_Node(4), 1);
s       = p(1);

% constraint and stationary residual
h = full(self.FuncObj.h(z));
c = full(self.FuncObj.c(z));
g = full(self.FuncObj.g(z, s));
LAG_grad = full(self.FuncObj.LAG_grad(Y, p));

% scaling parameter
scaling_dual = max([self.Option.KKT_scaling_max,...
    norm([gamma_h; gamma_c; gamma_g], 1)/(self.NLP.Dim.h + self.NLP.Dim.c + self.NLP.Dim.g)])/self.Option.KKT_scaling_max;
scaling_complementary = max([self.Option.KKT_scaling_max,...
    norm([gamma_c; gamma_g], 1)/(self.NLP.Dim.c + self.NLP.Dim.g)])/self.Option.KKT_scaling_max;

% KKT error
KKT_error_primal = norm([h;...
    min([zeros(self.NLP.Dim.c, 1), c], [], 2);...
    min([zeros(self.NLP.Dim.g, 1), g], [], 2)], inf);

KKT_error_dual = norm([LAG_grad';...
    min([zeros(self.NLP.Dim.c, 1), gamma_c], [], 2);...
    min([zeros(self.NLP.Dim.g, 1), gamma_g], [], 2)], inf);
KKT_error_dual = KKT_error_dual/scaling_dual;

KKT_error_complementary = norm([c .* gamma_c; g .* gamma_g], inf);
KKT_error_complementary = KKT_error_complementary/scaling_complementary;

KKT_error_total = max([KKT_error_primal, KKT_error_dual, KKT_error_complementary]);

end