function [tau_curr,rate_vec] = rate_opt_fun(opt_rates,opt_indices,K,c1,c0,activator_flag)   
    edge_key = [7 8 5 6 3 4 1 2];
    % establish determined edges    
    stable_indices = find(~ismember(1:8,opt_indices));
    % generate rate vector
    rate_vec = NaN(1,8);
    rate_vec(opt_indices) = opt_rates;
    if all(ismember(opt_indices,1:4))
        rate_vec(stable_indices) = rate_vec(edge_key(stable_indices));
    else
        rate_vec(stable_indices) = rate_vec(stable_indices-4) ...
            .* rate_vec(stable_indices-6) ./ rate_vec(stable_indices-2);  
    end
    % calculate mean
    m1 = fourStateProduction(rate_vec,c1,activator_flag);
    m0 = fourStateProduction(rate_vec,c0,activator_flag);
    % calculate noise
    s1 = sqrt(fourStateVariance(rate_vec,c1,activator_flag));
    s0 = sqrt(fourStateVariance(rate_vec,c0,activator_flag));      
    % calculate expected drift and diffusion rates
    V = .25*((m1-m0)^2*(s0^2+s1^2))/(s0^2*s1^2);
    D = .25*((m1-m0)^2*(s0^6+s1^6))/(s0^4*s1^4);
    % calcualte logL expectation    
    B = V*K/D;
    tau_curr = K/(2*V*sinh(B)) * (exp(B) + exp(-B) - 2);