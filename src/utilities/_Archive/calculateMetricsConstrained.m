function [metric_vec, profile, dd_prof, var_vec] = ...
                calculateMetricsConstrained(rate_vec,c_range,activator_flag)
            
%calculate expected derivative of pd rates wrt c
if activator_flag
    dmdcInv = calculateGradientActivator(rate_vec,1);
else
    dmdcInv = abs(calculateGradientRepressor(rate_vec,1));
end

% get c increment
dc = (c_range(2)-c_range(1))/2;

% c_range = c_range;
% calculate profile-wide metric
profile = fourStateProduction(rate_vec,c_range,activator_flag);
dd_prof = abs(diff(profile)/dc/2);
var_vec = fourStateVariance(rate_vec,c_range,activator_flag);

% calculate at center point
dmdc = abs(1/dmdcInv);
varPoint = fourStateVariance(rate_vec,1,activator_flag);

% binary metrics
ind1 = 101;
ind0 = 1;
c1 = c_range(ind1);
c0 = c_range(ind0);
pd1 = profile(ind1);
pd0 = profile(ind0);
fidelity = abs(log(pd1/pd0)/log(c0/c1));
dRange = abs((pd1-pd0)/(c1-c0));

% calculate max fidelity
p_vec = abs(log(profile(1) ./ profile));
c_vec = abs(log(c_range(1) ./ c_range));
f_vec = p_vec./c_vec;
max_i = randsample(1:numel(f_vec),1);

% calculate average decision time
m_critical = profile(51);
delta_m = abs(profile - m_critical);
decision_times = 2* var_vec ./ delta_m.^2;
t_critical = mean(decision_times([1:69]));

% calculate projected error rate at 40 minutes
error_vec = NaN(size(profile));
for p = 1:numel(profile)
    error_vec(p) = normcdf(-delta_m(p),0,sqrt(var_vec(p)/1200));
end

average_error = nanmean(error_vec([1:50]))+realmin;

metric_vec = [c_vec(max_i), p_vec(max_i), fidelity, dRange, dmdc, -.5*log(varPoint), -log(t_critical), log(average_error)];
