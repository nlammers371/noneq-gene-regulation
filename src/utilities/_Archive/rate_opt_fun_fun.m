function [tau_curr, rate_vec] = rate_opt_fun_fun(opt_rates,opt_indices,K,c1,c0,activator_flag)
    ref_indices = 1:8;
    stable_indices = ref_indices(~ismember(ref_indices,opt_indices));
    rate_vec = NaN(1,8);
    rate_vec(opt_indices) = opt_rates;
    rate_vec(stable_indices) = rate_vec(stable_indices-4);
    
%     % helper sub-functions
    m1_fun = @(rate_vec) fourStateProduction(rate_vec,c1,activator_flag);
    m0_fun = @(rate_vec) fourStateProduction(rate_vec,c0,activator_flag);

    s1_fun = @(rate_vec) sqrt(fourStateVariance(rate_vec,c1,activator_flag));
    s0_fun = @(rate_vec) sqrt(fourStateVariance(rate_vec,c0,activator_flag));

    V_fun = @(rate_vec) .25*((m1_fun(rate_vec)-m0_fun(rate_vec))^2*(s0_fun(rate_vec)^2+...
        s1_fun(rate_vec)^2))/(s0_fun(rate_vec)^2*s1_fun(rate_vec)^2);
    D_fun = @(rate_vec) .25*((m1_fun(rate_vec)-m0_fun(rate_vec))^2*(s0_fun(rate_vec)^6+s1_fun(rate_vec)^6))...
        /(s0_fun(rate_vec)^4*s1_fun(rate_vec)^4);


    % objective function
    tau_fun = @(rate_vec) K/(2*V_fun(rate_vec)*sinh(V_fun(rate_vec) * K / D_fun(rate_vec))) ...
        * (exp(V_fun(rate_vec) * K / D_fun(rate_vec)) + exp(-V_fun(rate_vec) * K / D_fun(rate_vec)) - 2);
    
    % calculate 
%     tau_curr = tau_fun(rate_vec);
%     m1 = fourStateProduction(rate_vec,c1,activator_flag);
%     m0 = fourStateProduction(rate_vec,c0,activator_flag);
% 
%     s1 = sqrt(fourStateVariance(rate_vec,c1,activator_flag));
%     s0 = sqrt(fourStateVariance(rate_vec,c0,activator_flag));
%   
%     %
%     % calculate expected drift and diffusion rates
%     V = .25*((m1-m0)^2*(s0^2+s1^2))/(s0^2*s1^2);
%     D = .25*((m1-m0)^2*(s0^6+s1^6))/(s0^4*s1^4);
%     % calcualte logL expectation    
%     B = V*K/D;
%     tau_curr = K/(2*V*sinh(B)) * (exp(B) + exp(-B) - 2);
    tau_curr = tau_fun(rate_vec);