function [metric_vec, rate_vec_new,metric_names,metric_ub_vec,metric_lb_vec] = calculateMetricsGeneral(rate_vec,h_max_flag)

if ~isempty(rate_vec)            
    % solve for c value where production rate is at half max
    rate_vec_new = rate_vec;
    c_val = 1;
    
    if h_max_flag
        c_val_init = fourStateHalfMaxGeneral(rate_vec_new);        
        rate_vec_new(2) = rate_vec_new(2)*c_val_init;
        rate_vec_new(6) = rate_vec_new(6)*c_val_init;    
    end
    
    % calculate production rate and sharpness     
    ProductionRate = fourStateProductionGeneral(rate_vec_new,c_val);    
    Sharpness = fourStateSharpnessGeneral(rate_vec_new,c_val); 
    
    % calculate thermodynamic force
    flux = log(prod(rate_vec_new(1:4))/prod(rate_vec_new(5:8)));
    
    % calculate variance
    VarPoint = fourStateVarianceGeneral(rate_vec_new,c_val);
    
    % calculate cycle time
    CycleTime = fourStateCycleTimeGeneral(rate_vec_new,c_val);            
    
    % generate vector to output
    metric_vec = [Sharpness, ProductionRate, flux, 1/sqrt(VarPoint), Sharpness/sqrt(VarPoint),CycleTime,c_val];
        
else
    metric_vec = [];  
    rate_vec_new = [];
end
metric_names = {'Sharpness','Production Rate','Flux','Precision','Information','CycleTime','c_val'};

% specify default bounds
metric_ub_vec = [Inf               Inf          Inf      Inf                 Inf                   Inf  Inf];
metric_lb_vec = [-Inf             -Inf         -Inf     -Inf                 -Inf                 -Inf -Inf];
