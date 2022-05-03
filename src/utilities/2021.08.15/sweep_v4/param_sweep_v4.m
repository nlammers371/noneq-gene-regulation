function simulation_results = param_sweep_v4(metric_indices,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Set parameter defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specify sampling hyperparameters
n_sim = 5; % number of independent runs to undertake
n_iters = 2e1; % iterations per sim
rate_bounds = log10([1e-4, 1e4]); % constrain transition rate magnitude
prop_sigma = 0.25;
n_seeds = 5; % number of new "seed" networks to generate per point on boundary
n_bound_samples = 20; % number of boundary samples to take
% n_b_points_max = 500; % maximum number of points to include in boundary
equilibrium_flag = 0; % if 1, detail balanced is enforced
rnd_seed = 'shuffle';

off_indices = [4 8]; % specify indives of unbinding rates
constrain_off_rates = true; % indicate whether to impose constrain on unbinding rates
off_val = 1;

% names of metric options
[~,metric_names,metric_ub_vec,metric_lb_vec] = calculateMetrics_v4([]);

% check for optional inputs
for i = 1:numel(varargin)
   if ischar(varargin{i}) && i <  numel(varargin)
       eval([varargin{i} ' = varargin{i+1};'])
   end
end

% calculate useful quantities
% n_samples_init = n_seeds*n_bound_samples;
rep_size = n_seeds*n_bound_samples;
  
tic
% iterate through specified model specs
simulation_results = struct;

% initialize arrays
rate_array = NaN(n_seeds*n_bound_samples*n_iters,8,n_sim);  
metric_array = NaN(n_seeds*n_bound_samples*n_iters,length(metric_names),n_sim);    

% record relevant hyperparameters        
simulation_results.rate_bounds = rate_bounds;     
simulation_results.metric_names = metric_names; 
simulation_results.edge_metric_indices = metric_indices;  
simulation_results.metric_lb_vec = metric_lb_vec;
simulation_results.metric_ub_vec = metric_ub_vec;

% initialize random number generator
rng(rnd_seed);

% define indexing vectors for convenience
rate_index = 1:8;

% iterate through specified number of independent runs   
for nti = 1:n_sim
  rates_temp = NaN(n_seeds*n_bound_samples*n_iters,8);
  metrics_temp = NaN(n_seeds*n_bound_samples*n_iters,length(metric_names));
  % generate initial array of samples
  lb_array = repmat(rate_bounds(1)/(5*prop_sigma),n_bound_samples*n_seeds,8);
  ub_array = repmat(rate_bounds(2)/(5*prop_sigma),n_bound_samples*n_seeds,8);
  new_rates = reshape(10.^(5*prop_sigma*trandn(lb_array,ub_array)),[],8); % function to sample from truncated normal dist
  if constrain_off_rates
    new_rates(:,off_indices) = off_val;
  end
  % enforce detailed balance if necessary
  if equilibrium_flag
    % first try just renormalizing backward rates
    net_flux_vec = prod(new_rates(:,1:4),2) ./ prod(new_rates(:,5:8),2);
    if ~constrain_off_rates
      new_rates(:,5:8) = net_flux_vec.^.25 .* new_rates(:,5:8);
    else
      rate_ft = ismember(rate_index,5:8)&~ismember(rate_index,off_indices);
      new_rates(:,rate_ft) = net_flux_vec.^(1/3) .* new_rates(:,rate_ft);
    end
    % find cases where this didn't work
    err_indices = find(max(log10(new_rates(:,5:8))<rate_bounds(1)|log10(new_rates(:,5:8))>rate_bounds(2),[],2));
    % resample these cases
    for e = err_indices'
      if ~constrain_off_rates
        k_sum = sum(log10(new_rates(e,1:4)));
        b_log = randfixedsum(4,1,k_sum,rate_bounds(1),rate_bounds(2));
        new_rates(e,5:8) = 10.^b_log;
      else
        k_sum = sum(log10(new_rates(e,1:4))) - log10(off_val);
        b_log = randfixedsum(3,1,k_sum,rate_bounds(1),rate_bounds(2));
        new_rates(e,rate_ft) = 10.^b_log;
      end
    end
  end
  % record
  rates_temp(1:rep_size,:) = new_rates;
  % calculate metrics  
  metrics_temp(1:rep_size,:) = calculateMetrics_v4(rates_temp(1:n_seeds*n_bound_samples,:));
  
  % now perform iterative edge-biased sampling
  for i = 2:n_iters     
    last_index = (i-1)*n_seeds*n_bound_samples;
    % extract and reshape arrays
    rate_array_curr = rates_temp(1:last_index,:);
    metric_array_curr = metrics_temp(1:last_index,metric_indices);    
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
    % generate expanded array for sampling
    orig_rate_array = log10(repmat(rate_array_curr(sample_indices,:),n_seeds,1));
    
    % generate variants
    lb_array = (rate_bounds(1)-orig_rate_array)/prop_sigma;
    ub_array = (rate_bounds(2)-orig_rate_array)/prop_sigma;
    new_rates = 10.^(orig_rate_array+reshape(prop_sigma*trandn(lb_array,ub_array),[],8));
    if constrain_off_rates
      new_rates(:,off_indices) = off_val;
    end
    % enforce detailed balance if necessary
    if equilibrium_flag        
      % first try just renormalizing backward rates
      net_flux_vec = prod(new_rates(:,1:4),2) ./ prod(new_rates(:,5:8),2);
      if ~constrain_off_rates
        new_rates(:,5:8) = net_flux_vec.^.25 .* new_rates(:,5:8);
      else
        rate_ft = ismember(rate_index,5:8)&~ismember(rate_index,off_indices);
        new_rates(:,rate_ft) = net_flux_vec.^(1/3) .* new_rates(:,rate_ft);
      end
      % find cases where this didn't work
      err_indices = find(max(log10(new_rates(:,5:8))<rate_bounds(1)|log10(new_rates(:,5:8))>rate_bounds(2),[],2));
      % resample these cases
      for e = err_indices'
        if ~constrain_off_rates
          k_sum = sum(log10(new_rates(e,1:4)));
          b_log = randfixedsum(4,1,k_sum,rate_bounds(1),rate_bounds(2));
          new_rates(e,5:8) = 10.^b_log;
        else
          k_sum = sum(log10(new_rates(e,1:4))) - log10(off_val);
          b_log = randfixedsum(3,1,k_sum,rate_bounds(1),rate_bounds(2));
          new_rates(e,rate_ft) = 10.^b_log;
        end
      end
    end
    
    % record
    rates_temp(last_index+1:i*n_seeds*n_bound_samples,:) = new_rates;
    % calculate metrics    
    metrics_temp(last_index+1:i*n_seeds*n_bound_samples,:) = calculateMetrics_v4(rates_temp(last_index+1:i*n_seeds*n_bound_samples,:));
  end  
  metric_array(:,:,nti) = metrics_temp;
  rate_array(:,:,nti) = rates_temp;
end       

% reshape output
simulation_results.metric_array = metric_array;        
simulation_results.rate_array = rate_array;    

toc
% save
% fpath = [outPath save_name];
% save(fpath,'simulation_results','-v7.3')