% script to use edge sampling to explore performance landscape of
% proofreading networks
clear
close all
% specify paths
addpath('../utilities')
outPath = '../../out/parameter_sweep_general/';
mkdir(outPath)

% define specific edges 
m_space = [1 7; 3 5];
% specify sampling hyperparameters
n_mcmc = 5; % number of independent runs to undertake
n_iters = 5e2;
metric_indices = [1 2];
n_bound_points = 150;
c_range = linspace(2/11,20/11,101);
% set save name
save_name = ['edge_sim_results_N' sprintf('%04d',n_iters)  '.mat'];

% sampling params and hyperparams
rate_bounds = [1e-2, 1e2]; % constrain transition rate magnitude
n_seeds = 10; % number of "seed" spots to generate per point on boundary
% iterate through specified model specs
simulation_results = struct;
for mdi = 1:size(m_space,1)  
    activator_flag = mdi==1;
    % determine edge groups
    binding_edges = m_space(mdi,:);    
    % initialize arrays
    rate_array = NaN(1+n_seeds*n_iters,8,n_mcmc);        
    % metric array     
    metric_array_eq = NaN(1+n_seeds*n_iters,3,n_mcmc);    
    metric_array_noneq = NaN(1+n_seeds*n_iters,3,n_mcmc);    
    % record relevant hyperparameters        
    simulation_results(mdi).model_spec = binding_edges;                
    % iterate through specified number of independent MCMC runs
    tic
    for mci = 1:n_mcmc                         
        % initialize rates 
        rate_vec_init = randsample(logspace(log10(rate_bounds(1)),log10(rate_bounds(2))),8,true);           
        % decide on initial "move"
        orig_rate_vec = sample_rates_general(rate_vec_init, rate_bounds);
        rate_array(1,:,mci) = orig_rate_vec;                 
        % calculate non-eq metrics
        metric_vec = calculateMetricsConstrained(orig_rate_vec,c_range,activator_flag);        
        metric_array_noneq(1,:,mci) = metric_vec;                                    
        % eq rates
        eq_rate_vec = orig_rate_vec;
        eq_rate_vec(5:8) = randsample(eq_rate_vec(1:4),4,false);
        metric_vec_eq = calculateMetricsConstrained(eq_rate_vec,c_range,activator_flag);        
        metric_array_eq(1,:,mci) = metric_vec_eq;                                    
        % generate first round of mutants            
        for k = 1:n_seeds
            % sample eq rates
            new_rate_vec = sample_rates_general(orig_rate_vec, rate_bounds);
            % record
            rate_array(1+k,:,mci) = new_rate_vec;
            % metrics
            metric_vec = calculateMetricsConstrained(new_rate_vec,c_range, activator_flag);                              
            metric_array_noneq(1+k,:,mci) = metric_vec;
            % impose detailed balance
            eq_rate_vec = new_rate_vec;
            eq_rate_vec(5:8) = randsample(eq_rate_vec(1:4),4,false);
            metric_vec_eq = calculateMetricsConstrained(eq_rate_vec,c_range,activator_flag);        
            metric_array_eq(1+k,:,mci) = metric_vec_eq;                                    
        end                
        % set initial specs for edge sampling
        bound_points = 1:k+1;
        prev_point = metric_array_noneq(1,metric_indices,mci);
        % perform edge sampling
        for n = 2:n_iters                  
            % find set of boundary points from which to take sample          
            use_indices = bound_points(~isnan(metric_array_noneq(bound_points,1,mci)));
            n_sample = min([numel(use_indices),n_bound_points]);
            use_indices = randsample(use_indices,n_sample,false);
            % filter metric array
            metric_array_curr = metric_array_noneq(use_indices,metric_indices,mci);
            metric_array_curr = normalize([metric_array_curr(:,1) log(metric_array_curr(:,2))]);
            % calculate boundary points
            bound_points = unique(use_indices(boundary(metric_array_curr,1)));  
            % randomly select next point...up-weight points that are
            % distant from current point
            diff_array = metric_array_noneq(bound_points,metric_indices,mci) - prev_point;
            dist_array = sqrt(sum((diff_array ./ mean(diff_array)).^2,2));            
            mt_index = randsample(bound_points,1,true,dist_array);    
            % extract metrics and corresponding rates
            prev_point = metric_array_noneq(mt_index,metric_indices,mci);                       
            rates_orig = rate_array(mt_index,:,mci);            
            % generate mutants            
            for k = 1:n_seeds
                ind = 1+(n-1)*n_seeds+k;
                % sample rates
                new_rate_vec = sample_rates_general(rates_orig, rate_bounds);
                rate_array(ind,:,mci) = new_rate_vec;               
                % calculate eq performance metrics
                metric_vec = calculateMetricsConstrained(rate_array(ind,:,mci),c_range,activator_flag);                                   
                metric_array_noneq(ind,:,mci) = metric_vec;                                           
                % impose detailed balance
                eq_rate_vec = new_rate_vec;
                eq_rate_vec(5:8) = randsample(eq_rate_vec(1:4),4,false);
                metric_vec_eq = calculateMetricsConstrained(eq_rate_vec,c_range,activator_flag);        
                metric_array_eq(ind,:,mci) = metric_vec_eq;                                    
                % add points
                bound_points = [bound_points ind];
            end                              
        end
        disp(['completed simulation ' sprintf('%03d',mci) ' of ' num2str(n_mcmc)...
            ' (model ' sprintf('%02d',mdi) ' of ' num2str(size(m_space,1)) ')'])        
    end
    toc
    simulation_results(mdi).metric_array_eq = metric_array_eq;
    simulation_results(mdi).metric_array_noneq = metric_array_noneq;        
    simulation_results(mdi).rate_array = rate_array;
end

save([outPath save_name],'simulation_results','-v7.3')