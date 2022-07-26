function f = computeStateEquation(plant, tau, x, p)
%computeStateEquation
%   compute State Equation

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

