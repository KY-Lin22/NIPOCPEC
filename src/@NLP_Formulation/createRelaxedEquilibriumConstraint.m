function [constraint_without_param, constraint_with_param] = createRelaxedEquilibriumConstraint(self, OCPEC)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

import casadi.*
% variable
lambda = SX.sym('lambda', OCPEC.Dim.lambda, 1);
eta = SX.sym('eta', OCPEC.Dim.lambda, 1);
% parameter
s = SX.sym('s', 1, 1);

% formulate constraint
switch OCPEC.VISetType
    case 'box_constraint'
        constraint_without_param_formula = [...
            lambda - OCPEC.bl;...
            OCPEC.bu - lambda];
        constraint_with_param_formula = [...
            s  - (lambda - OCPEC.bl) .* eta;...
            s  + (OCPEC.bu - lambda) .* eta];

        constraint_without_param = Function('constraint_without_param',...
            {lambda, eta}, {constraint_without_param_formula},...
            {'lambda', 'eta'}, {'constraint_without_param'});
        constraint_with_param = Function('constraint_with_param',...
            {lambda, eta, s}, {constraint_with_param_formula},...
            {'lambda', 'eta', 's'}, {'constraint_with_param'});

    case 'nonnegative_orthant'
        constraint_without_param_formula = [lambda; eta];
        constraint_with_param_formula = s - lambda .* eta;

        constraint_without_param = Function('constraint_without_param',...
            {lambda, eta}, {constraint_without_param_formula},...
            {'lambda', 'eta'}, {'constraint_without_param'});
        constraint_with_param = Function('constraint_with_param',...
            {lambda, eta, s}, {constraint_with_param_formula},...
            {'lambda', 'eta', 's'}, {'constraint_with_param'});

end

end