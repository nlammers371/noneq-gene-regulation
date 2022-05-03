function [metric_vec, metric_names] = ...
                calculateMetricsV2(rate_vec,c_range,frequency_flag)
            
if ~isempty(rate_vec)            
%     % calculate expected derivative of pd rates wrt c
%     if frequency_flag
%         dmdcInv = calculateGradientActivator(rate_vec,1);
%     else
%         dmdcInv = abs(calculateGradientRepressor(rate_vec,1));
%     end

    % calculate at center point
%     dMdC = abs(1/dmdcInv);
    VarPoint = fourStateVariance(rate_vec,1,frequency_flag);
    % calculate profile-wide metric
    dc = c_range(2)-c_range(1);
    profile = fourStateProduction(rate_vec,c_range,frequency_flag);
    dd_prof = abs(diff(profile)/dc/2);
    dMdC = dd_prof(numel(profile)/2);
    ProductionRate = mean(profile([numel(profile)/2 numel(profile)/2+1]));
    % fold metrics
    ind1 = numel(c_range);
    ind0 = 1;
    c1 = c_range(ind1);
    c0 = c_range(ind0);
    pd1 = profile(ind1);
    pd0 = profile(ind0);
    fidelity = abs(log(pd1/pd0)/log(c0/c1));
    dRange = abs((pd1-pd0)/(c1-c0));
    
    % randomly sample fidelity 
    p_vec = abs(log(profile(1) ./ profile));
    c_vec = abs(log(c_range(1) ./ c_range));
    f_vec = p_vec./c_vec;
    max_i = randsample(1:numel(f_vec),1);
       
    % randomly sample fidelity    
    dp_vec = profile(end) - profile(1:end-1);
    dc_vec = c_range(end) - c_range(1:end-1);
    s_vec = dp_vec./dc_vec;
    index_s = randsample(1:numel(s_vec),1);
        
    % calculate average decision time
    % m_critical = profile(51);
    % delta_m = abs(profile - m_critical);
    % decision_times = 2* var_vec ./ delta_m.^2;
    % t_critical = mean(decision_times([1:69]));
    % 
    % % calculate projected error rate at 40 minutes
    % error_vec = NaN(size(profile));
    % for p = 1:numel(profile)
    %     error_vec(p) = normcdf(-delta_m(p),0,sqrt(var_vec(p)/1200));
    % end
    % 
    % average_error = nanmean(error_vec([1:50]))+realmin;
    metric_vec = [c_vec(max_i), p_vec(max_i),dc_vec(index_s), dp_vec(index_s), fidelity, dRange, dMdC, sqrt(VarPoint),... 
        ProductionRate,dMdC/sqrt(VarPoint)];
else
    metric_vec = [];
end
metric_names = {'fold-concentration','fold-production','diff-concentration','diff-production',...
    'fidelity','DiscreteSharpness','Sharpness','Noise','Production Rate','Positional Error'};
