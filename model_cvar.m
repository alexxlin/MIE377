function  cvar = model_cvar(weights, scenarios, C_I)
    
    %Compute VaR using CI 
    percentile = size(scenarios,1) * C_I;
    percentile = int16(percentile);
    
    %actual loss for each scenario based on optimal weights
    f_xy = weights' * scenarios' * -1;
    
    %find the average loss of the part of loss that exceeds VaR
    sorted_fxy = sort(f_xy);
    cvar = mean(sorted_fxy(percentile+1:end));
end