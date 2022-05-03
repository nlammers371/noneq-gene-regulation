function metric_vec = calculatePerformanceMetrics(rate_vec,spec_edges,c)

k_vec0 = rate_vec(1:4);
r_vec0 = rate_vec(5:8);
% make high vectors
k_vec1 = k_vec0;
k_vec1(spec_edges(1)) = k_vec1(spec_edges(1))*c;
r_vec1 = r_vec0;
r_vec1(spec_edges(2)) = r_vec1(spec_edges(2)-4)*c;
% calculate expected production rates 
pd0 = fourStateProduction(k_vec0,r_vec0);
pd1 = fourStateProduction(k_vec1,r_vec1);
fidelity = abs(log10(pd1/pd0) / log10(c));
% calculate asymptotic variance
asv0 = fourStateVariance(k_vec0,r_vec0);
asv1 = fourStateVariance(k_vec1,r_vec1);
% now calculate expected logL divergence rate
information_rate = log10(exp(1))*.25*((pd1-pd0)/sqrt(asv1))^2 + .25*((pd1-pd0)/sqrt(asv0))^2;
% calculate inverse of expecte time for high and low states to diverge
c_time = ((pd1-pd0)/(sqrt(asv1) + sqrt(asv0)))^-2;
% production rate
mean_rate = mean([pd0 pd1]);
% information rate
% information_rate = .5*((asv0*(pd1-pd0)^2)/(2*asv1) + (asv1*(pd1-pd0)^2)/(2*asv0) - 1);
% concatenate
metric_vec = [fidelity, information_rate, c_time, sign(pd1-pd0), mean_rate];