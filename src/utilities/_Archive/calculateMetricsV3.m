function [metric_vec, metric_names] = calculateMetricsV3(rate_vec,repressor_flag,h_max_flag)

if ~isempty(rate_vec)            
    % solve for c value where production rate is at half max
    if h_max_flag
        c_val = fourStateHalfMax(rate_vec,repressor_flag);
    else
        c_val = 1;
    end
    % calculate production rate and sharpness     
    ProductionRate = fourStateProduction_v2(rate_vec,c_val,repressor_flag);    
    Sharpness = fourStateSharpness_v2(rate_vec,c_val,repressor_flag); 
    % aclculate thermodynamic force
    flux = log(prod(rate_vec(1:4))/prod(rate_vec(5:8)));
    % calculate variance
    VarPoint = fourStateVariance_v2(rate_vec,c_val,repressor_flag);
    TauOff = fourStateTauOff(rate_vec,c_val,repressor_flag);    
    Pi2 = fourStatePi2(rate_vec,c_val,repressor_flag);    
    metric_vec = [Sharpness*c_val, ProductionRate, flux, ProductionRate./sqrt(VarPoint), abs(Sharpness)/sqrt(VarPoint),TauOff,Pi2];
else
    metric_vec = [];
end
metric_names = {'Sharpness','Production Rate','Flux','Variance','Information','ReturnTime','Pi2'};
