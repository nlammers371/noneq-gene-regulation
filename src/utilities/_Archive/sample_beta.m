function [rates_noneq, beta_val_new] = sample_beta(rate_bounds, rate_vec_orig, beta_edge, beta_orig)
    % sample beta                             
    beta_bounds = [1, rate_bounds(2)/rate_vec_orig(beta_edge)];
    pd = makedist('Lognormal',log(beta_orig),2);
    pd_truncated = truncate(pd,beta_bounds(1),beta_bounds(2));
    beta_val_new = random(pd_truncated,1);
    % calculate noneq metrics and record
    rates_noneq = rate_vec_orig;                
    rates_noneq(beta_edge) = rates_noneq(beta_edge)*beta_val_new;