function K = computeVI_Function(plant, tau, x, p)
%computeVI_Function
%   compute function K of VI

K = autoGen_K(tau, x, p);

end

