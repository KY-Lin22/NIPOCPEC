function f = computeStateEquation_Function(plant, tau, x, p)
%computeStateEquation_Function
%   compute State Equation Function

if plant.computeInvM
    q_Dim = round(1/2 * plant.Dim.x);
    M = autoGen_M(x);
    H = autoGen_H(tau, x, p);
    f = [x(q_Dim + 1 : end);...
         M\H];   
else
    f = autoGen_f(tau, x, p);   
end

end

