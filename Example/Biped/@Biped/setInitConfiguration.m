function q = setInitConfiguration(plant, x, xita_torso, xita_thigh_1, xita_calf_1, xita_thigh_2)
%setInitConfiguration
%   DoubleSupport
Dim_q = 1/2 * plant.Dim.x;
q = zeros(Dim_q, 1);

q(1) = x;
q(3) = xita_torso;
q(4) = xita_thigh_1;
q(5) = xita_calf_1;
q(6) = xita_thigh_2;

q(2) = plant.linkLength(2) * cos(q(4)) + plant.linkLength(3) * cos(q(5)); % z

keel_to_ground = q(2) - plant.linkLength(2) * cos(q(6));
q(7) = - acos(keel_to_ground/ plant.linkLength(3)); % leg 2 (ee of calf 2)
end

