% script to use edge sampling to explore performance landscape of
% proofreading networks
clear
close all
% specify paths
addpath('../utilities')
outPath = '../../out/parameter_sweep_constrained/';
mkdir(outPath)

% define specific edges 
m_space = [1 7; 3 5];
% specify sampling hyperparameters
n_mcmc = 1; % number of independent runs to undertake
n_iters = 1e3;
metric_indices = [5 6];
n_bound_points = 500;
c_range = linspace(2/11,20/11,101);
n_samples_init = 500;
full_mod = 10;
rate_bounds = [1e-4, 1e2]; % constrain transition rate magnitude
n_seeds = 5; % number of "seed" spots to generate per point on boundary
% set save name
metric_names = {'pointwise_error','average_error','fidelity','info_rate','dmdc','InvSigma','fidelity_orig','info_rate_orig'};
% set save name
save_name = ['edge_sim_results_' metric_names{metric_indices(1)} '_' ...
    metric_names{metric_indices(2)}  '.mat'];

% iterate through specified model specs
simulation_results = struct;
for mdi = 1:size(m_space,1)  
    activator_flag = mdi==1;
    % determine edge groups
    binding_edges = m_space(mdi,:);    
    % initialize arrays
    rate_array = NaN(n_samples_init+n_seeds*n_iters,8,8,n_mcmc);    
    beta_edge_array = repmat(1:8,n_samples_init+n_seeds*n_iters,1,n_mcmc);
    beta_factor_array = NaN(n_samples_init+n_seeds*n_iters,8,n_mcmc);
    % metric array 
    metric_array_eq = NaN(n_samples_init+n_seeds*n_iters,numel(metric_names),8,n_mcmc);    
    metric_array_noneq = NaN(n_samples_init+n_seeds*n_iters,numel(metric_names),8,n_mcmc);    
    % record relevant hyperparameters        
    simulation_results(mdi).model_spec = binding_edges;                
    simulation_results(mdi).N = n_iters*n_seeds + n_samples_init;                
    simulation_results(mdi).rate_bounds = rate_bounds;
    % iterate through specified number of independent MCMC runs
    tic
    for mci = 1:n_mcmc     
        for beta_edge = 1:8            
            for n = 1:n_samples_init
                % initialize rates 
                rate_vec_init = randsample(logspace(log10(rate_bounds(1)),log10(rate_bounds(2))),4,true);           
                % decide on initial "move"
                new_rate_vec = sample_rates_constrained(rate_vec_init, rate_bounds);
                rate_array(n,:,beta_edge,mci) = new_rate_vec;                 
                % calculate non-eq metrics
                metric_vec_eq = calculateMetricsConstrained(new_rate_vec,c_range,activator_flag);        
                metric_array_eq(n,:,beta_edge,mci) = metric_vec_eq;  
                % sample Beta                  
                beta_val = max(1,rand()*rate_bounds(2)/new_rate_vec(beta_edge));       
                new_rate_vec_beta = new_rate_vec;
                new_rate_vec_beta(beta_edge) = new_rate_vec_beta(beta_edge)*beta_val;
                % record            
                beta_factor_array(n,beta_edge,mci) = beta_val;
                % calculate non-eq metrics
                metric_vec_neq = calculateMetricsConstrained(new_rate_vec_beta,c_range,activator_flag);  
                metric_array_noneq(n,:,beta_edge,mci) = metric_vec_neq;                
            end            
            tic            
            % set initial specs for edge sampling
            bound_points = 1:n_samples_init;            
            % perform edge sampling
            for n = 1:n_iters                                  
                % filter metric array
                if mod(n,full_mod) == 0
                    use_indices = find(~isnan(metric_array_noneq(:,1,beta_edge,mci)))';                                        
                else                    
                    n_sample = min([numel(bound_points),n_bound_points]);
                    use_indices = randsample(bound_points,n_sample,false);                    
                end            
                metric_array_curr = metric_array_noneq(use_indices,metric_indices,beta_edge,mci);
                metric_array_curr = normalize(metric_array_curr);
                % calculate boundary points                 
                bpoints = unique(boundary(metric_array_curr,1));
                bound_points = use_indices(bpoints);  
                curr_points = metric_array_curr(bpoints,:);
                D = pdist2(curr_points,curr_points);
                D(eye(size(D,1))==1) = Inf;
                dist_vec = min(D,[],2);
                % randomly select next point...up-weight points that are in sparse regions                      
                mt_index = randsample(bound_points,1,true,dist_vec);  
                % extract metrics and corresponding rates                
                rates_orig = rate_array(mt_index,:,beta_edge,mci);
                beta_orig = beta_factor_array(mt_index,beta_edge,mci);                
                % generate mutants            
                for k = 1:n_seeds
                    ind = n_samples_init+(n-1)*n_seeds+k;
                    % sample rates
                    seed_rate_vec = sample_rates_constrained(rates_orig, rate_bounds);                    
                    rate_array(ind,:,beta_edge,mci) = seed_rate_vec;               
                    % calculate eq performance metrics
                    metric_vec_eq = calculateMetricsConstrained(rate_array(ind,:,beta_edge,mci),c_range,activator_flag);                                   
                    metric_array_eq(ind,:,beta_edge,mci) = metric_vec_eq;                           
                    % sample beta
                    [rates_noneq, beta_val_new] = sample_beta(rate_bounds, seed_rate_vec, beta_edge, beta_orig);
                    beta_factor_array(ind,beta_edge,mci) = beta_val_new;
                    % calculate non-eq performance metrics
                    metric_vec_noneq = calculateMetricsConstrained(rates_noneq, c_range,activator_flag);                                   
                    metric_array_noneq(ind,:,beta_edge,mci) = metric_vec_noneq; 
                    % add to point list         
                    bound_points = [bound_points ind];
                end                  
            end
            toc
        end
        disp(['completed simulation ' sprintf('%03d',mci) ' of ' num2str(n_mcmc)...
            ' (model ' sprintf('%02d',mdi) ' of ' num2str(size(m_space,1)) ')'])        
    end
    toc
    simulation_results(mdi).metric_array_eq = metric_array_eq;
    simulation_results(mdi).metric_array_noneq = metric_array_noneq;
    simulation_results(mdi).beta_edge_array = beta_edge_array;
    simulation_results(mdi).beta_factor_array = beta_factor_array;
    simulation_results(mdi).rate_array = rate_array;
end

save([outPath save_name],'simulation_results','-v7.3')