function setEqualityConstraints(OCPEC, C)
%setEqualityConstraints
%   Detailed explanation goes here

if ~isempty(C)
    OCPEC.C = [OCPEC.C;...
               C];
    OCPEC.Dim.eta = size(OCPEC.C, 1);
    OCPEC.eta = sym('eta', [OCPEC.Dim.eta, 1]);
end

end

