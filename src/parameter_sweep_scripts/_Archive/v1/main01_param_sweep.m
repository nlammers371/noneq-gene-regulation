% script to use MCMC sampling to explore performance landscape of
% proofreading networks
clear
close all
% specify paths
addpath('../utilities')
outPath = '../../out/parameter_search/';
mkdir(outPath)

% define specific edges 
m_space = [1 3; 3 1];
% specify sampling hyperparameters
model_spec_vec = [1:2]; 
n_mcmc = 1; % number of independent runs to undertake
n_iters = 1e4;
equilibrium_flag = 0; % constrained to equilibrium?
metric_names = {'cooperativity', 'information_rate', 'critical time', 'production_sign', 'production_rate'};
metric_indices = [1 2];
n_bound_points = 150;
% set save name
save_name = ['edge_sim_results_eq' num2str(equilibrium_flag) '_' metric_names{metric_indices(1)} '_' metric_names{metric_indices(2)} '.mat'];

% sampling params and hyperparams
prior_weights = ones(1,8)/8; % controls initial size of each edg
rate_bounds = [1e-2, 1e2]; % constrain transition rate magnitude
delta_mu_in = 10; % ratio of high-to-low concentrations (assumed fixed)

% initialize beta distribution to sample rate value multipliers
proposal_vec = logspace(-2,2,1000);
norm_edge = 8;
% number of "seed" spots to generate per point on boundary
n_seeds = 10;
% iterate through specified model specs
simulation_results = struct;
for mdi = 1:numel(model_spec_vec)
    % set eq norm edge
    norm_edge_iter = norm_edge;
    % determine specific edges
    spec_edges = m_space(model_spec_vec(mdi),:);
    unbinding_edges = [spec_edges(2) spec_edges(1)+4];
    if norm_edge_iter == unbinding_edges(2)
        options = 5:8;
        options = options(options~=unbinding_edges(2));
        norm_edge_iter = options(1);
    end
    % initialize arrays
    rate_array = NaN(1+n_seeds*n_iters,8,n_mcmc);
    % metric array 
    metric_array = NaN(1+n_seeds*n_iters,numel(metric_names),n_mcmc);    
    % record relevant hyperparameters    
    simulation_results(mdi).model_id = model_spec_vec(mdi);    
    simulation_results(mdi).model_spec = spec_edges;        
    simulation_results(mdi).proposal_vec = proposal_vec;
    simulation_results(mdi).equilibrium_flag = equilibrium_flag;
    simulation_results(mdi).metric_names = metric_names;
    simulation_results(mdi).metric_indices = metric_indices;
    % iterate through specified number of independent MCMC runs
    tic
    for mci = 1:n_mcmc
        r_indices = 5:8;
        % initialize rates 
        rate_vec_init = rand(1,8)*5;   
        rate_vec_init(unbinding_edges) = 1;
        if equilibrium_flag
            rate_vec_init(norm_edge_iter) = prod(rate_vec_init(1:4)) / prod(rate_vec_init(r_indices(r_indices~=norm_edge_iter)));  
        end                   
        rate_array(1,:,mci) = rate_vec_init;        
        % calculate performance_metrics
        metric_vec = calculatePerformanceMetrics(delta_mu_in,spec_edges,...
            rate_vec_init(1:4),rate_vec_init(5:8));        
        metric_array(1,:,mci) = metric_vec;
        % generate first round of mutants
        for k = 1:n_seeds
            prop_factors = randsample(proposal_vec,8,true);
            prop_factors(unbinding_edges) = 1;
            prop_rates = rate_vec_init.*prop_factors;
            % impose boundary conditions and equilibrium conditions as
            % needed
            boundary_violations = prop_rates>rate_bounds(2) | prop_rates<rate_bounds(1);
            prop_rates(boundary_violations) = rate_vec_init(boundary_violations);
            if equilibrium_flag
                prop_rates(norm_edge_iter) = prod(prop_rates(1:4)) / prod(prop_rates(r_indices(r_indices~=norm_edge_iter)));
            end    
            % second check
            boundary_violations = prop_rates>rate_bounds(2) | prop_rates<rate_bounds(1);                
            if any(boundary_violations) && equilibrium_flag
                continue
            elseif any(boundary_violations)
                error('wtf')
            end
            rate_array(1+k,:,mci) = prop_rates;
            metric_vec = calculatePerformanceMetrics(delta_mu_in,spec_edges, ...
                rate_array(1+k,1:4), rate_array(1+k,5:8));      
            metric_array(1+k,:,mci) = metric_vec;
        end       
        bound_points = 1:k+1;
        prev_point = metric_array(1,metric_indices);
        % perform MCMC sampling
        for n = 2:n_iters               
            % find set of boundary points from which to take sample          
            use_indices = bound_points(~isnan(metric_array(bound_points,1)));
            n_sample = min([numel(use_indices),n_bound_points]);
            use_indices = randsample(use_indices,n_sample,false);
            % filter metric array
            metric_array_curr = metric_array(use_indices,metric_indices,mci);
            metric_array_curr = metric_array_curr ./ nanmean(metric_array_curr);            
            bound_points = unique(use_indices(boundary(metric_array_curr,1)));   
            diff_array = metric_array(bound_points,metric_indices) - prev_point;
            dist_array = sqrt(sum((diff_array ./ mean(diff_array)).^2,2));            
            mt_index = randsample(bound_points,1,true,dist_array);    
            prev_point = metric_array(mt_index,metric_indices);
            % randomly pick one of the outlier points            
            rates_orig = rate_array(mt_index,:,mci);
            % generate mutants            
            for k = 1:n_seeds
                ind = 1+(n-1)*n_seeds+k;
                prop_factors = randsample(proposal_vec,8,true);
                prop_factors(unbinding_edges) = 1;
                prop_rates = rates_orig.*prop_factors;
                % impose boundary conditions and equilibrium conditions as
                % needed                                
                boundary_violations = prop_rates>rate_bounds(2) | prop_rates<rate_bounds(1);
                prop_rates(boundary_violations) = rates_orig(boundary_violations);
                if equilibrium_flag
                    prop_rates(norm_edge_iter) = prod(prop_rates(1:4)) / prod(prop_rates(r_indices(r_indices~=norm_edge_iter)));
                end    
                % second check
                boundary_violations = prop_rates>rate_bounds(2) | prop_rates<rate_bounds(1);                
                if any(boundary_violations) && equilibrium_flag
                    continue
                elseif any(boundary_violations)
                    error('wtf')
                end
                rate_array(ind,:,mci) = prop_rates;
                metric_vec = calculatePerformanceMetrics(delta_mu_in,spec_edges,...
                    rate_array(ind,1:4), rate_array(ind,5:8));                                   
                metric_array(ind,:,mci) = metric_vec;                
                bound_points = [bound_points ind];
            end            
        end             
        disp(['completed simulation ' sprintf('%03d',mci) ' of ' num2str(n_mcmc)...
            ' (model ' sprintf('%02d',mdi) ' of ' num2str(numel(model_spec_vec)) ')'])        
    end
    toc
    simulation_results(mdi).metric_array = metric_array;
    simulation_results(mdi).rate_array = rate_array;
end

save([outPath save_name],'simulation_results','-v7.3')