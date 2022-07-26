function setDynVarLimit(plant, tau_Max, tau_Min, x_Max, x_Min)
%setDynVarLimit
%   Detailed explanation goes here
if ~isempty(tau_Max)
    if all(size(tau_Max) == [plant.Dim.tau, 1])
        plant.DynVarLimit.tau_Max = tau_Max;
    else
        error('please specify tau_Max with the correct dim')
    end
end

if ~isempty(tau_Min)
    if all(size(tau_Min) == [plant.Dim.tau, 1])
        plant.DynVarLimit.tau_Min = tau_Min;
    else
        error('please specify tau_Min with the correct dim')
    end
end

if ~isempty(x_Max)
    if all(size(x_Max) == [plant.Dim.x, 1])
        plant.DynVarLimit.x_Max = x_Max;
    else
        error('please specify x_Max with the correct dim')
    end
end

if ~isempty(x_Min)
    if all(size(x_Min) == [plant.Dim.x, 1])
        plant.DynVarLimit.x_Min = x_Min;
    else
        error('please specify x_Min with the correct dim')
    end
end

end

