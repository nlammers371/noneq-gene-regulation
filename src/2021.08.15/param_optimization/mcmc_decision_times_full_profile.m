clear 
close all
addpath('../utilities')
outPath = '../../out/mcmc_optimization/';
mkdir(outPath);
figPath = '../../fig/mcmc_optimization/';
mkdir(figPath);

% generate input profile data
x_vec = linspace(0,1,25);
f_diff = 1.1;
scale_factor = log(f_diff) / (x_vec(2)-x_vec(1));
c_prefactor = exp(scale_factor/2);
c_profile = c_prefactor*exp(-x_vec.*scale_factor);

% iterate through neighboring points and calculate expected decision times
% under different constraints

% params
rate_key = [1 2 3 4; 
            7 8 5 6];
edge_opt_cell = {4,5,6,7,8,5:8};
activator_flag = 0;
K = log(100);
rate_bounds = [1e-4,1e1];
prop_sigma = .2;
n_samples = 500;
ref_indices = 1:8;
temp = 100;
burn_in = 200;
% initialize arrays
decision_time_array = NaN(numel(c_profile)-1,n_samples,numel(edge_opt_cell));
opt_rate_array = NaN(n_samples,8,numel(edge_opt_cell));
for r = 1:numel(edge_opt_cell)
    opt_indices = unique([1:4 edge_opt_cell{r}]);
    nr = numel(opt_indices);
    stable_indices = find(~ismember(rate_key(2,:),opt_indices));            
    % set optimization routine params
    options = logspace(-2,1);
    init_rates = randsample(options,nr);
    % define caller function
%     rate_fun_call = @(rate_vec) rate_opt_fun(rate_vec,opt_indices,K,c1,c0,activator_flag);
    curr_tau_vec = NaN(1,numel(c_profile)-1);
    for c = 1:numel(c_profile)-1
        c1 = c_profile(c);
        c0 = c_profile(c+1);
        [curr_tau_vec(c), current_rate_vec] = rate_opt_fun(init_rates,opt_indices,K,c1,c0,activator_flag);
    end
    % perform MCMC sampling        
    for n = 1:n_samples   
        % record current values        
        decision_time_array(:,n,r) = curr_tau_vec;
        mean_tau_curr = nanmean(curr_tau_vec);
        opt_rate_array(n,:,r) = current_rate_vec;        
        % sample new rates
        new_rates = sample_rates_general(current_rate_vec(opt_indices),rate_bounds,prop_sigma);  
        new_rate_vec = NaN(1,8);
        new_rate_vec(opt_indices) = new_rates;
        new_rate_vec(rate_key(2,stable_indices)) = new_rates(rate_key(1,stable_indices));        
        % calculate new decision time
        new_tau_vec = NaN(1,numel(c_profile)-1);
        parfor c = 1:numel(c_profile)-1
            c1 = c_profile(c);
            c0 = c_profile(c+1);
            [new_tau_vec(c), ~] = rate_opt_fun(new_rate_vec(opt_indices),opt_indices,K,c1,c0,activator_flag);
        end
        mean_tau_new = nanmean(new_tau_vec);
        % MH move
        score = exp((mean_tau_curr-mean_tau_new)/temp);
        update = score > rand();
        if update
            current_rate_vec = new_rate_vec;
            curr_tau_vec = new_tau_vec;
        end            
    end
end

sim_struct = struct;
sim_struct.edge_opt_vec = edge_opt_cell;
sim_struct.c_profile = c_profile;
sim_struct.rate_bounds = rate_bounds;
im_struct.activator_flag = activator_flag;
sim_struct.logL_rate = K;
sim_struct.decision_time_array = decision_time_array;
sim_struct.opt_rate_array = opt_rate_array;
% save 
save([outPath 'rate_optimization_full_profile_act' num2str(activator_flag) '.mat'],'sim_struct');