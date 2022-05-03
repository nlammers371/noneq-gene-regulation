% script to examine how limits of system behavior change as a function of
% energy input
clear
close all
addpath('../utilities')

%% (1) Information Rate
c_val_info = 1;
ub_vec = 1e4*ones(size(rates_init));
lb_vec = 1e-4*ones(size(rates_init));
options = optimoptions('fmincon','Display','off');

% first solve for overall max value. This will serve as a scale parameter
% for constraint component
rates_init = rand(1,8);
info_fun = @(rate_vec) -abs(fourStateSharpnessGeneral(rate_vec,c_val_info))/sqrt(fourStateVarianceGeneral(rate_vec,c_val_info));
flux_fun = @(rate_vec) abs(log(prod(rate_vec(1:4))/prod(rate_vec(5:8))));
info_rates = fmincon(info_fun,rates_init,[],[],[],[],lb_vec,ub_vec,[],options);
max_info = -ob_fun(info_rates);
max_flux = flux_fun(info_rates);

%% calculate twmperature vector
T1 = 1e-2;
T2 = 1e2*max_flux/max_info;

temperature_vec = logspace(log10(T1),log10(T2),50);
flux_vec = NaN(1,numel(temperature_vec));
info_vec = NaN(1,numel(temperature_vec));
rate_array = NaN(numel(temperature_vec),8);

for t = 1:numel(temperature_vec)
    tic
    rates_init = rand(1,8);
    % define weighted objective function 
    info_ob_fun = @(rate_vec) info_fun(rate_vec) + flux_fun(rate_vec)/temperature_vec(t);
    info_rates = fmincon(info_ob_fun,rates_init,[],[],[],[],lb_vec,ub_vec,[],options);
    info_vec(t) = -info_fun(info_rates);
    flux_vec(t) = flux_fun(info_rates);
    rate_array(t,:) = info_rates;
    toc
end




  
         