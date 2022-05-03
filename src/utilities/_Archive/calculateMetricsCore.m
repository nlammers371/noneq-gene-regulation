function metric_vec = ...
                calculateMetricsCore(rate_vec,c_vec,activator_flag)
            
%calculate expected derivative of pd rates wrt c
if activator_flag
    dmdcInv = calculateGradientActivator(rate_vec,1);
else
    dmdcInv = abs(calculateGradientRepressor(rate_vec,1));
end

% calculate at center point
dmdc = abs(1/dmdcInv);
varPoint = fourStateVariance(rate_vec,1,activator_flag);

% binary metrics
c1 = c_vec(2);
c0 = c_vec(1);
pd1 = fourStateProduction(rate_vec,c_vec(1),activator_flag);
pd0 = fourStateProduction(rate_vec,c_vec(2),activator_flag);
fidelity = abs(log(pd1/pd0)/log(c0/c1));
dRange = abs((pd1-pd0)/(c1-c0));


metric_vec = [fidelity, dRange, dmdc/sqrt(varPoint)];
