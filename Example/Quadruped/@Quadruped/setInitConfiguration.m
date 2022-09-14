function q = setInitConfiguration(plant, x, xita)
%setInitConfiguration
%   Detailed explanation goes here

Dim_q = 1/2 * plant.Dim.x;
q = zeros(Dim_q, 1);

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

