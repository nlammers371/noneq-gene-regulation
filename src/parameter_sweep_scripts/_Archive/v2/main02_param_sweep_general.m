% script to use edge sampling to explore performance landscape of
% proofreading networks
clear
close all
% specify paths
addpath('../utilities')
outPath = '../../out/parameter_search_v2/';
mkdir(outPath)

% define specific edges 
m_space = [1 7; 3 5];
% specify sampling hyperparameters
n_mcmc = 1; % number of independent runs to undertake
n_iters = 5e2;
metric_indices = [10,11];
c_range = linspace(2/11,20/11,101);
rate_bounds = [1e-4, 1e2]; % constrain transition rate magnitude
proposal_vec = logspace(-1,1,1000);
n_b_points_max = 1000;
n_seeds = 5; % number of "seed" spots to generate per point on boundary
n_bound_samples = 20;
n_samples_init = n_seeds*n_bound_samples;
% names of metric options
metric_names = {'pointwise-error','average-error','dynamicRange','fidelity','info-rate','dmdc',...
    'InvSigma','decisionRate','boundPosition','fold-concentration','fold-production'};
% set save name
save_name = ['edge_sim_results_' metric_names{metric_indices(1)} '_' metric_names{metric_indices(2)} '.mat'];

% iterate through specified model specs
simulation_results = struct;
for mdi = 1:size(m_space,1)  
    activator_flag = mdi==1;
    % determine edge groups
    binding_edges = m_space(mdi,:);    
    % initialize arrays
    rate_array_noneq = NaN(n_seeds,8,n_bound_samples,n_iters,n_mcmc);        
    rate_array_eq = NaN(n_seeds,8,n_bound_samples,n_iters,n_mcmc);        
    % metric array     
    metric_array_eq = NaN(n_seeds,numel(metric_names),n_bound_samples,n_iters,n_mcmc);        
    metric_array_noneq = NaN(n_seeds,numel(metric_names),n_bound_samples,n_iters,n_mcmc);        
    % record relevant hyperparameters        
    simulation_results(mdi).model_spec = binding_edges;                    
    simulation_results(mdi).c_range = c_range;                
    % iterate through specified number of independent MCMC runs
    tic
    for mci = 1:n_mcmc         
        tic
        parfor n = 1:n_bound_samples
            for s = 1:n_seeds
                % decide on initial "move"            
                orig_rate_vec = randsample(logspace(log10(rate_bounds(1)),log10(rate_bounds(2))),8,true);                       
                rate_array_noneq(s,:,n,1,mci) = orig_rate_vec;                 
                % calculate non-eq metrics
                metric_vec = calculateMetricsConstrained(orig_rate_vec,c_range,activator_flag);        
                metric_array_noneq(s,:,n,1,mci) = metric_vec;              
                % eq rates
                eq_rate_vec = orig_rate_vec;
                eq_rate_vec(5:8) = randsample(eq_rate_vec(1:4),4,false);
                metric_vec_eq = calculateMetricsConstrained(eq_rate_vec,c_range,activator_flag);        
                metric_array_eq(s,:,n,1,mci) = metric_vec_eq;                 
                rate_array_eq(s,:,n,1,mci) = eq_rate_vec;
            end
        end
        toc       
        % perform edge sampling
        for n = 2:n_iters                                                                                                                               
            rate_array_curr = reshape(permute(rate_array_noneq(:,:,:,1:n-1,mci),[1 4 3 2]),[],8);
            metric_array_curr = reshape(permute(metric_array_noneq(:,metric_indices,:,1:n-1,mci),[1 4 3 2]),[],2);           
            metric_array_curr = normalize(metric_array_curr);      
            metric_array_curr = metric_array_curr - min(metric_array_curr);
            index_vec = 1:numel(metric_array_curr(:,1));
            % estimate point densities
            [N,Xedges,Yedges,binX,binY] = histcounts2(metric_array_curr(:,1),metric_array_curr(:,2));
            % calculate center of mass
%             edge_wt_array = 1./bwdist(N==0);
            lin_indices = sub2ind(size(N),binX,binY);
            density_weights = N(lin_indices).^-2;
%             edge_weights = edge_wt_array(lin_indices);
            if mod(n,1) == 0
                sample_indices = randsample(index_vec,n_bound_samples,true,density_weights);  
            else
                b_samp_points = randsample(index_vec,min(n_b_points_max,numel(index_vec)),true,density_weights.^-2);
                bound_indices = b_samp_points(boundary(metric_array_curr(b_samp_points,1),metric_array_curr(b_samp_points,2),1));
                sample_indices = randsample(bound_indices, n_bound_samples, true);
            end
          
            parfor b = 1:n_bound_samples
                mt_index = sample_indices(b);
                % extract metrics and corresponding rates           
                rates_orig = rate_array_curr(mt_index,:);            
                % generate mutants            
                for k = 1:n_seeds                    
                    % sample rates
                    prop_factors = randsample(proposal_vec,8,true);
                    new_rate_vec = prop_factors.*rates_orig;
                    new_rate_vec(new_rate_vec>rate_bounds(2)) = rate_bounds(2);
                    new_rate_vec(new_rate_vec<rate_bounds(1)) = rate_bounds(1);
                    rate_array_noneq(k,:,b,n,mci) = new_rate_vec;               
                    % calculate eq performance metrics
                    metric_vec = calculateMetricsConstrained(new_rate_vec,c_range,activator_flag);                                   
                    metric_array_noneq(k,:,b,n,mci) = metric_vec;                                           
                    % impose detailed balance
                    eq_rate_vec = new_rate_vec;                    
                    eq_rate_vec(5:8) = randsample(eq_rate_vec(1:4),4,false);                    
                    metric_vec_eq = calculateMetricsConstrained(eq_rate_vec,c_range,activator_flag);        
                    metric_array_eq(k,:,b,n,mci) = metric_vec_eq;          
                    rate_array_eq(k,:,b,n,mci) = eq_rate_vec;                                                                 
                end  
            end            
        end
        disp(['completed simulation ' sprintf('%03d',mci) ' of ' num2str(n_mcmc)...
            ' (model ' sprintf('%02d',mdi) ' of ' num2str(size(m_space,1)) ')'])        
    end
    toc
    simulation_results(mdi).metric_array_eq = metric_array_eq;
    simulation_results(mdi).metric_array_noneq = metric_array_noneq;        
    simulation_results(mdi).rate_array_noneq = rate_array_noneq;
    simulation_results(mdi).rate_array_eq = rate_array_eq;
end

save([outPath save_name],'simulation_results','-v7.3')