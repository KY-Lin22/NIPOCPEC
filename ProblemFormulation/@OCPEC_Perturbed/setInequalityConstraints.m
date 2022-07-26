function setInequalityConstraints(OCPEC, G)
%setInequalityConstraints
%   Detailed explanation goes here

if ~isempty(G)    
    OCPEC.G = [OCPEC.G;...
              G];
    OCPEC.Dim.sigma = size(OCPEC.G, 1);
    OCPEC.sigma = sym('sigma', [OCPEC.Dim.sigma, 1]);
end

end

