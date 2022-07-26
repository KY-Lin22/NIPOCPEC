function K = computeVIFunc(plant, tau, x, p)
%computeVIFunc
%   compute function K of VI

K = autoGen_K(tau, x, p);

end

