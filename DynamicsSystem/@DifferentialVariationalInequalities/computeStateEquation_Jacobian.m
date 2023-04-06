function [ftau, fx, fp] = computeStateEquation_Jacobian(plant, tau, x, p)
%computeStateEquation_Jacobian
%   compute State Equation Jacobian

Dim = plant.Dim;
          
switch plant.computeStateEquationMethod
    case 1
        [ftau, fx, fp] = autoGen_ftau_fx_fp(tau, x, p);
    case 2
        q_Dim = round(1/2 * Dim.x);
        M = autoGen_M(x);
        H = autoGen_H(tau, x, p);
        Mx = autoGen_Mx(x);
        [Htau, Hx, Hp] = autoGen_Htau_Hx_Hp(tau, x, p);
        %
        InvMH = M\H;
        dInvMH_dtau = M\Htau;
        MxInvMH = zeros(q_Dim, Dim.x);
        for i = 1 : Dim.x
            Mx_i = Mx(:, (i - 1)*q_Dim + 1 : i * q_Dim);
            MxInvMH(:, i) = Mx_i * InvMH;
        end
        dInvMH_dx = M\(Hx - MxInvMH);
        dInvMH_dp = M\Hp;
        
        ftau = [zeros(q_Dim, Dim.tau);...
            dInvMH_dtau];
        fx = [zeros(q_Dim, q_Dim), eye(q_Dim);...
            dInvMH_dx];
        fp = [zeros(q_Dim, Dim.p);...
            dInvMH_dp];
        
    case 3
        q_Dim = plant.qDim;
        M = autoGen_M(x);
        H = autoGen_H(tau, x, p);
        Mx = autoGen_Mx(x);
        [Htau, Hx, Hp] = autoGen_Htau_Hx_Hp(tau, x, p);   
        [f_giventau, f_givenx, f_givenp] = autoGen_ftau_fx_fp(tau, x, p);
        %
        InvMH = M\H;
        dInvMH_dtau = M\Htau;
        MxInvMH = zeros(q_Dim, Dim.x);
        for i = 1 : Dim.x
            Mx_i = Mx(:, (i - 1)*q_Dim + 1 : i * q_Dim);
            MxInvMH(:, i) = Mx_i * InvMH;
        end
        dInvMH_dx = M\(Hx - MxInvMH);
        dInvMH_dp = M\Hp;
        
        ftau = [zeros(q_Dim, Dim.tau);...
            dInvMH_dtau;...
            f_giventau];
        fx = [zeros(q_Dim, q_Dim), eye(q_Dim), zeros(q_Dim, Dim.x - 2 * q_Dim);...
            dInvMH_dx;...
            f_givenx];
        fp = [zeros(q_Dim, Dim.p);...
            dInvMH_dp;...
            f_givenp];        
end

end

