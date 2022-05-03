% script to use MCMC sampling to explore optimal network topologies
clear
close all
% specify paths
addpath('../utilities')
outPath = '../../out/parameter_search/';
mkdir(outPath)

% define specific edges 
m_space = [1 3; 3 1];
% specify sampling hyperparameters
model_spec_vec = [1:2]; % 5 corresponds to Hopfield/Ninio architecture
n_mcmc = 10; % number of independent runs to undertake
n_iters = 1e4;
equilibrium_flag = 1; % constrained to equilibrium?
metric_names = {'cooperativity', 'information_rate', 'critical_time', 'production_sign', 'production_rate'};
metric_index = 1;
proposal_vec = logspace(-.1,.1,1000);
% set save name
save_name = ['mcmc_opt_results_eq' num2str(equilibrium_flag) '_' metric_names{metric_index} '.mat'];

% sampling params and hyperparams
rate_bounds = [1e-2, 1e4]; % constrain transition rate magnitude
r_abs = rate_bounds(2) / rate_bounds(1);
delta_mu_in = 10; % ratio of high-to-low concentrations (assumed fixed)
n_mvrs = log(r_abs);
n_factor = rate_bounds(end);
beta_vec = [70 10];
beta = beta_vec(metric_index); % controls "temperature"
% iterate through specified model specs
mcmc_results = struct;
for mdi = 1:numel(model_spec_vec)    
    % determine specific edges
    spec_edges = m_space(model_spec_vec(mdi),:);
    unbinding_edges = [spec_edges(2) spec_edges(1)+4];   
    k_sample_indices = find((1:4)~=unbinding_edges(1));
    r_sample_indices = 4 + find((5:8)~=unbinding_edges(2));
    % initialize arrays
    rate_array = NaN(n_iters,8,n_mcmc);    
    % metric array 
    metric_array = NaN(n_iters,numel(metric_names),n_mcmc);     
    % record relevant hyperparameters    
    mcmc_results(mdi).model_id = model_spec_vec(mdi);    
    mcmc_results(mdi).model_spec = spec_edges;            
    mcmc_results(mdi).equilibrium_flag = equilibrium_flag;
    mcmc_results(mdi).metric_names = metric_names;
    mcmc_results(mdi).metric_index = metric_index;
    % iterate through specified number of independent MCMC runs
    tic
    for mci = 1:n_mcmc      
        % sample multiplicative factors that set scale for k and r cycles
        prop_factors = randsample(logspace(log10(rate_bounds(1)+100),log10(rate_bounds(2)-100),1000),8);
        prop_factors(unbinding_edges) = 1;
        if equilibrium_flag            
            prop_factors(r_sample_indices) = randsample(prop_factors(k_sample_indices),3,false);
        end
        rate_vec_init = prop_factors;                
        rate_array(1,:,mci) = rate_vec_init;        
        % calculate performance_metrics
        metric_vec = calculatePerformanceMetrics...
            (delta_mu_in,spec_edges,rate_vec_init(1:4),rate_vec_init(5:8));        
        metric_array(1,:,mci) = metric_vec;
        % perform MCMC sampling
        ct = 0;
        for n = 2:n_iters               
            % find set of boundary points from which to take sample                                      
            rates_orig = rate_array(n-1,:,mci);                        
            metrics_orig = metric_array(n-1,:,mci);
            % generate proposed move
            prop_factors = randsample(proposal_vec,8);            
            if equilibrium_flag                
                prop_factors(r_sample_indices) = randsample(prop_factors(k_sample_indices),3,false);
            end   
            rate_vec_prop = prop_factors .* rates_orig;
            rate_vec_prop(unbinding_edges) = 1;
            bound_flag = any(rate_vec_prop>rate_bounds(2) | rate_vec_prop<rate_bounds(1));
            if ~equilibrium_flag
                rate_vec_prop(rate_vec_prop>rate_bounds(2)) = rate_bounds(2);
                rate_vec_prop(rate_vec_prop<rate_bounds(1)) = rate_bounds(1);                                             
            end
            % calculate metric values
            metric_vec_prop = calculatePerformanceMetrics...
            (delta_mu_in,spec_edges,rate_vec_prop(1:4),rate_vec_prop(5:8));                        
            if equilibrium_flag && bound_flag
                metric_vec_prop(metric_index) = -Inf;
            end
            % compare scores and decide
            [~, mi] = max([rand() exp(beta*(metric_vec_prop(metric_index)-metrics_orig(metric_index)))]);
            if mi == 1
                rate_array(n,:,mci) = rates_orig;                
                metric_array(n,:,mci) = metrics_orig;
                ct = ct + 1;
            elseif mi == 2
                rate_array(n,:,mci) = rate_vec_prop;                
                metric_array(n,:,mci) = metric_vec_prop;
                ct = 0;
            else
                error('uh oh')
            end
            if ct > 500
                error('hmmmm')
            end
        end             
        disp(['completed simulation ' sprintf('%03d',mci) ' of ' num2str(n_mcmc)...
            ' (model ' sprintf('%02d',mdi) ' of ' num2str(numel(model_spec_vec)) ')'])        
    end
    toc
    % generate longform arrays
    metric_array_perm = permute(metric_array,[1 3 2]);
    metric_array_long = reshape(metric_array_perm,[],size(metric_array,2));
    rate_array_perm = permute(rate_array,[1 3 2]);
    rate_array_long = reshape(rate_array_perm,[],size(rate_array,2));
    
    metric_table = array2table(metric_array_long,'VariableNames',metric_names);
    mcmc_results(mdi).metric_array = metric_table;
    mcmc_results(mdi).rate_array = rate_array_long;
end

save([outPath save_name],'mcmc_results','-v7.3')