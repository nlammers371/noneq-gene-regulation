% script to use edge sampling to explore performance landscape of
% proofreading networks
clear
close all
% specify paths
addpath('../utilities')
outPath = '../../out/parameter_sweep_analyses/';
mkdir(outPath)

% define specific edges 
m_space = [1 7; 3 5];
% specify sampling hyperparameters
n_mcmc = 5; % number of independent runs to undertake
n_iters = 1e3;
metric_indices = [1 3];
n_bound_points = 500;
c_range = linspace(2/11,20/11,101);
decision_points = [1 2];
metric_names = {'pointwise_error','average_error','binary_fidelity','binary_info_rate','dmdc','varPoint'};
% set save name
save_name = ['edge_sim_results_eq' metric_names{metric_indices(1)} '_' ...
    metric_names{metric_indices(2)} '_N' sprintf('%04d',n_iters)  '.mat'];

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
    metric_array_eq = NaN(1+n_seeds*n_iters,numel(metric_names),n_mcmc);    
    % record relevant hyperparameters        
    simulation_results(mdi).model_spec = binding_edges;                
    % iterate through specified number of independent MCMC runs
    tic
    for mci = 1:n_mcmc                         
        % initialize rates 
        new_rate_vec = randsample(logspace(log10(rate_bounds(1)),log10(rate_bounds(2))),8,true);           
        % decide on initial "move"
%         new_rate_vec = sample_rates_general(rate_vec_init, rate_bounds);
        % impose detailed balance
        balance_factor = prod(new_rate_vec(1:4))/prod(new_rate_vec(5:8));
        prop_vec = new_rate_vec(5:8)*balance_factor;
        options = find(prop_vec>=rate_bounds(1)&prop_vec<=rate_bounds(2));
        if ~isempty(options)
            if numel(options) == 1
                norm_edge = options+4;
            else
                norm_edge = randsample(4+options,1);
            end
            new_rate_vec(norm_edge) = new_rate_vec(norm_edge)*balance_factor;
        else
            new_rate_vec(5:8) = randsample(new_rate_vec(1:4),4,false);
        end        
        rate_array(1,:,mci) = new_rate_vec;                 
        % calculate non-eq metrics
        metric_vec = calculateMetricsConstrained(new_rate_vec,c_range,activator_flag);        
        metric_array_eq(1,:,mci) = metric_vec;                                                  
        % generate first round of mutants            
        for k = 1:n_seeds
            % sample eq rates
            seed_rate_vec = sample_rates_general(new_rate_vec, rate_bounds);
            % impose detailed balance
            balance_factor = prod(seed_rate_vec(1:4))/prod(seed_rate_vec(5:8));
            prop_vec = seed_rate_vec(5:8)*balance_factor;
            options = find(prop_vec>=rate_bounds(1)&prop_vec<=rate_bounds(2));
            if ~isempty(options)
                if numel(options) == 1
                    norm_edge = options+4;
                else
                    norm_edge = randsample(4+options,1);
                end                
                seed_rate_vec(norm_edge) = seed_rate_vec(norm_edge)*balance_factor;
            else
                seed_rate_vec(5:8) = randsample(seed_rate_vec(1:4),4,false);
            end  
            if any(seed_rate_vec>rate_bounds(2)|seed_rate_vec<rate_bounds(1))
                error('problem with bounds')
            end
            % record
            rate_array(1+k,:,mci) = seed_rate_vec;
            % metrics
            metric_vec = calculateMetricsConstrained(seed_rate_vec,c_range, activator_flag);                              
            metric_array_eq(1+k,:,mci) = metric_vec;                                               
        end                
        % set initial specs for edge sampling
        bound_points = 1:k+1;
        prev_point = metric_array_eq(1,metric_indices,mci);
        % perform edge sampling
        for n = 2:n_iters                  
            % find set of boundary points from which to take sample          
            use_indices = bound_points(~isnan(metric_array_eq(bound_points,1,mci)));
            n_sample = min([numel(use_indices),n_bound_points]);
            use_indices = randsample(use_indices,n_sample,false);
            % filter metric array
            metric_array_curr = metric_array_eq(use_indices,metric_indices,mci);
            metric_array_curr = normalize([metric_array_curr(:,1) log(metric_array_curr(:,2))]);                                    
            
            % calculate boundary points
            bpoints = unique(boundary(metric_array_curr,1));
            bound_points = use_indices(bpoints);  
            curr_points = metric_array_curr(bpoints,:);
            D = pdist2(curr_points,curr_points);
            D(eye(size(D,1))==1) = Inf;
            dist_vec = min(D,[],2);
            % randomly select next point...up-weight points that are in sparse regions                      
            mt_index = randsample(bound_points,1,true,dist_vec.^2); 
                        
            rates_orig = rate_array(mt_index,:,mci);            
            % generate mutants            
            for k = 1:n_seeds
                ind = 1+(n-1)*n_seeds+k;
                % sample rates
                new_rate_vec = sample_rates_general(rates_orig, rate_bounds);
                % impose detailed balance
                balance_factor = prod(new_rate_vec(1:4))/prod(new_rate_vec(5:8));
                prop_vec = new_rate_vec(5:8)*balance_factor;
                options = find(prop_vec>=rate_bounds(1)&prop_vec<=rate_bounds(2));
                if ~isempty(options)
                    if numel(options) == 1
                        norm_edge = options+4;
                    else
                        norm_edge = randsample(4+options,1);
                    end
                    new_rate_vec(norm_edge) = new_rate_vec(norm_edge)*balance_factor;
                else
                    new_rate_vec(5:8) = randsample(new_rate_vec(1:4),4,false);
                end  
                if any(new_rate_vec>rate_bounds(2)|new_rate_vec<rate_bounds(1))
                    error('problem with bounds')
                end
                rate_array(ind,:,mci) = new_rate_vec;               
                % calculate eq performance metrics
                metric_vec = calculateMetricsConstrained(rate_array(ind,:,mci),c_range,activator_flag);                                   
                metric_array_eq(ind,:,mci) = metric_vec;                                                                                          
                % add points
                bound_points = [bound_points ind];
            end                              
        end
        disp(['completed simulation ' sprintf('%03d',mci) ' of ' num2str(n_mcmc)...
            ' (model ' sprintf('%02d',mdi) ' of ' num2str(size(m_space,1)) ')'])        
    end
    toc
    simulation_results(mdi).metric_array_eq = metric_array_eq;    
    simulation_results(mdi).rate_array = rate_array;
end

save([outPath save_name],'simulation_results','-v7.3')