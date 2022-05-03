function [simulation_results, fpath] = param_sweep_v3(metric_indices,varargin)
addpath('../utilities')
outPath = '../../out/bivariate_parameter_sweeps_v3/';
mkdir(outPath)

%%%% Set parameter defaults
% specify sampling hyperparameters
n_sim = 5; % number of independent runs to undertake
n_iters = 2e1; % iterations per sim
rate_bounds = [1e-2, 1e2]; % constrain transition rate magnitude
proposal_vec = logspace(-.5,.5,1000);
n_seeds = 5; % number of "seed" spots to generate per point on boundary
n_bound_samples = 20;
n_b_points_max = 500;
h_max_flag = 0;
equilibrium_flag = 0;
% names of metric options
[~, ~,metric_names,metric_ub_vec,metric_lb_vec] = calculateMetricsGeneral([],[]);
% metric bounds
% check for optional inputs
for i = 1:numel(varargin)
   if ischar(varargin{i}) && i <  numel(varargin)
       eval([varargin{i} ' = varargin{i+1};'])
   end
end

% calculate useful quantities
n_samples_init = n_seeds*n_bound_samples;
mem_steps = ceil(n_b_points_max/n_samples_init);

% set save name
save_name = ['param_sweep_results_' metric_names{metric_indices(1)} '_' ...
    metric_names{metric_indices(2)} '_eq' num2str(equilibrium_flag) '_hm' num2str(h_max_flag) '.mat'];
  
tic
% iterate through specified model specs
simulation_results = struct;


% initialize arrays
rate_array = NaN(n_seeds,8,n_bound_samples,n_iters*n_sim);  
metric_array = NaN(n_seeds,numel(metric_names),n_bound_samples,n_iters*n_sim);    

% record relevant hyperparameters        
simulation_results.rate_bounds = rate_bounds;     
simulation_results.metric_names = metric_names; 
simulation_results.edge_metric_indices = metric_indices;  
simulation_results.metric_lb_vec = metric_lb_vec;
simulation_results.metric_ub_vec = metric_ub_vec;
simulation_results.h_max_flag = h_max_flag;

% iterate through specified number of independent MCMC runs   
restart_indices = 1:n_iters:n_sim*n_iters;
for nti = 1:n_iters*n_sim
    % if at a restart index, draw set of initial points
    if ismember(nti,restart_indices)
        for b = 1:n_bound_samples
            for s = 1:n_seeds
                pass_flag = false;
                while ~pass_flag
                    % decide on initial "move"                             
                    orig_rate_vec_raw = randsample(logspace(log10(rate_bounds(1)),log10(rate_bounds(2))),8,true);                                                             
                    % apply detailed balance if necessary
                    if equilibrium_flag 
                        flux = prod(orig_rate_vec_raw(1:4))/prod(orig_rate_vec_raw(5:8));
                        orig_rate_vec_raw(5:8) = flux^.25 * orig_rate_vec_raw(5:8);
                    end
                    % calculate  metrics
                    [metric_vec,orig_rate_vec] = calculateMetricsGeneral(orig_rate_vec_raw,h_max_flag);                     
                    pass_flag = all(metric_vec>=metric_lb_vec)&&all(metric_vec<=metric_ub_vec) ...
                         &&all(orig_rate_vec<=rate_bounds(2))&&all(orig_rate_vec>=rate_bounds(1));
                end
                rate_array(s,:,b,nti) = orig_rate_vec;                                                                                                    
                metric_array(s,:,b,nti) = metric_vec;                    
            end
        end  
    else                 
        % calculate lookback distance
        mem_max = mod(nti-1,n_iters);
        mem = min(mem_max,mem_steps);

        % extract and reshape arrays
        rate_array_curr = reshape(permute(rate_array(:,:,:,nti-mem:nti-1),[1 4 3 2]),[],8);
        metric_array_curr = reshape(permute(metric_array(:,metric_indices,:,nti-mem:nti-1),[1 4 3 2]),[],2);    
        use_indices = find(max(isnan(metric_array_curr),[],2)==0);

        % obtain normalized array containing current metric values                
        metric_array_norm = normalize(metric_array_curr(use_indices,:));   

        % find boundaries of metric space
        bpoints = boundary(metric_array_norm,.5);
        dist_array = pdist2(metric_array_norm(bpoints,:),metric_array_norm(bpoints,:));
        dist_array(eye(numel(bpoints))==1) = Inf;
        dist_vec = nanmin(dist_array,[],2);

        % resample boundary points according to nearest neighbor
        % distance        
        sample_indices = randsample(use_indices(bpoints),n_bound_samples,true,dist_vec);                      

        for b = 1:n_bound_samples
            mt_index = sample_indices(b);
            % extract metrics and corresponding rates           
            rates_orig = rate_array_curr(mt_index,:);                   
            % generate mutants            
            for k = 1:n_seeds          
                pass_flag = false;
                while ~pass_flag
                    % sample rates      
                    prop_factors = randsample(proposal_vec,8,true);                        
                    new_rate_vec_raw = prop_factors.*rates_orig;  
                    new_rate_vec_raw(new_rate_vec_raw>rate_bounds(2)) = rate_bounds(2);
                    new_rate_vec_raw(new_rate_vec_raw<rate_bounds(1)) = rate_bounds(1);     
                    % apply detailed balance if necessary
                    if equilibrium_flag 
                        flux = prod(new_rate_vec_raw(1:4))/prod(new_rate_vec_raw(5:8));
                        new_rate_vec_raw(5:8) = flux^.25 * new_rate_vec_raw(5:8);
                    end
                    % calculate eq performance metrics
                    [metric_vec, new_rate_vec] = calculateMetricsGeneral(new_rate_vec_raw,h_max_flag);       
%                     pass_flag = all(metric_vec>=metric_lb_vec)&all(metric_vec<=metric_ub_vec);%| ...
                        %any(isnan(metric_vec));
                    pass_flag = all(metric_vec>=metric_lb_vec)&&all(metric_vec<=metric_ub_vec) ...
                        && all(new_rate_vec<=rate_bounds(2)) && all(new_rate_vec>=rate_bounds(1));
                end
                % record                   
                rate_array(k,:,b,nti) = new_rate_vec;                                                                           
                metric_array(k,:,b,nti) = metric_vec;                                                                                                                            
            end  
        end  
    end        
end       
% reshape output
rate_array_out = reshape(permute(rate_array(:,:,:,:),[1 4 3 2]),[],8);
metric_array_out = reshape(permute(metric_array(:,:,:,:),[1 4 3 2]),[],numel(metric_names)); 
simulation_results.metric_array = metric_array_out;        
simulation_results.rate_array = rate_array_out;    

toc
% save
fpath = [outPath save_name];
save(fpath,'simulation_results','-v7.3')