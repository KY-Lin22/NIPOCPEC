function var_Sequence = TrajectoryInterpolation(var_Init, var_End, nStages)

% checking dim
if ~all(size(var_Init) == size(var_End))
    warning('wrong in setting var_Init and var_End');
end
% interpolation
var_Dim = size(var_Init, 1);
var_Sequence = zeros(var_Dim, nStages);
interpolationStepSize = 1/(nStages - 1);
for i = 1 : var_Dim
    var_Sequence(i, :) = interp1([var_Init(i), var_End(i)], 1:interpolationStepSize:2, 'linear');
end

end

