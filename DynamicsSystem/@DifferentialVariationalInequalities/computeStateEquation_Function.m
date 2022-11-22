function f = computeStateEquation_Function(plant, tau, x, p)
%computeStateEquation_Function
%   compute State Equation Function

switch plant.computeStateEquationMethod
    case 1
        f = autoGen_f(tau, x, p); 
        
    case 2
        q_Dim = round(1/2 * plant.Dim.x);
        M = autoGen_M(x);
        H = autoGen_H(tau, x, p);
        f = [x(q_Dim + 1 : q_Dim * 2);...
            M\H];
        
    case 3
        q_Dim = plant.qDim;
        M = autoGen_M(x);
        H = autoGen_H(tau, x, p);  
        f_given = autoGen_f(tau, x, p);
        f = [x(q_Dim + 1 : q_Dim * 2);...
            M\H;...
            f_given];
end

end

