function q = setInitConfiguration(plant, x, xita)
%setInitConfiguration
%   Detailed explanation goes here

q_Dim = plant.qDim;
q = zeros(q_Dim, 1);

q(1) = x;
q(3) = 1/2 * pi;
q(4) = -xita;
q(5) = xita;
q(6) = -xita;
q(7) = xita;
q(8) = -xita;
q(9) = xita;
q(10) = -xita;
q(11) = xita;
q(2) = plant.linkLength(2) * cos(q(4)) + plant.linkLength(3) * cos(q(5));

end

