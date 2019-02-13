function  cvar = cvar_value(weights, scenarios, C_I)
    threshold = int16(size(scenarios,1) * C_I);
    f_xy = -weights' * scenarios';
    sorted_fxy = sort(f_xy);
    var = sorted_fxy(threshold);
    cvar = mean(sorted_fxy(threshold+1:end));
end