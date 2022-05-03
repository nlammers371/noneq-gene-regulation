function [simulation_results, fpath] = param_sweep_fun(metric_indices,equilibrium_flag,varargin)
addpath('../utilities')
outPath = '../../out/bivariate_parameter_sweeps/';
mkdir(outPath)

%%%% Set parameter defaults

% specify sampling hyperparameters
n_sim = 5; % number of independent runs to undertake
n_iters = 2e1; % iterations per sim
c_range = linspace(2/2.1,2/2.1*1.1,100);
rate_bounds = [1e-4, 1e4]; % constrain transition rate magnitude
proposal_vec = logspace(-.5,.5,1000);
n_seeds = 5; % number of "seed" spots to generate per point on boundary
n_bound_samples = 20;
n_b_points_max = 500;

% check for optional inputs
for i = 1:numel(varargin)
   if ischar(varargin{i}) && i <  numel(varargin)
       eval([varargin{i} ' = varargin{i+1};'])
   end
end

% calculate useful quantities
n_samples_init = n_seeds*n_bound_samples;
mem_steps = ceil(n_b_points_max/n_samples_init);

% names of metric options
[~, metric_names] = calculateMetricsV2([],[],[]);
% set save name
save_name = ['param_sweep_results_' metric_names{metric_indices(1)} '_' ...
    metric_names{metric_indices(2)} 'eq' num2str(equilibrium_flag) '.mat'];
  
tic
% iterate through specified model specs
simulation_results = struct;

% initialize arrays
rate_array = NaN(n_seeds,8,n_bound_samples,n_iters*n_sim);  
metric_array = NaN(n_seeds,numel(metric_names),n_bound_samples,n_iters*n_sim); 

% record relevant hyperparameters                          
simulation_results(mdi).c_range = c_range; 
simulation_results(mdi).rate_bounds = rate_bounds;     
simulation_results(mdi).metric_names = metric_names; 
simulation_results(mdi).edge_metric_indices = metric_indices;     
 
% iterate through specified number of independent MCMC runs   

  restart_indices = 1:n_iters:n_sim*n_iters;
  for nti = 1:n_iters*n_sim
      % if at a restart index, draw set of initial points
      if ismember(nti,restart_indices)
          for b = 1:n_bound_samples
              for s = 1:n_seeds
                  % decide on initial "move"            
                  % impose detailed balance if necessary
                  if equilibrium_flag
                      exit_flag = false;
                      while ~exit_flag
                          orig_rate_vec = randsample(logspace(log10(rate_bounds(1)),log10(rate_bounds(2))),8,true);                                              
                          flux = prod(orig_rate_vec(1:4))/prod(orig_rate_vec(5:8));
                          orig_rate_vec(5:8) = flux^.25 * orig_rate_vec(5:8);
                          exit_flag = all(orig_rate_vec<=rate_bounds(2))&&all(orig_rate_vec>=rate_bounds(1));
                      end                            
                  else
                      orig_rate_vec = randsample(logspace(log10(rate_bounds(1)),log10(rate_bounds(2))),8,true);                                              
                  end
                  rate_array(s,:,b,nti) = orig_rate_vec;                                                                
                  % calculate  metrics
                  metric_vec = calculateMetricsV2(orig_rate_vec,c_range,frequency_flag);                                        
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
          use_indices = 1:size(metric_array_curr,1);
          if all(ismember(metric_indices,[4,6]))
              use_indices = find(metric_array_curr(:,metric_indices==4)>.05);
          end
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
              metric_vec_orig = calculateMetricsV2(rates_orig,c_range,frequency_flag);
              orig_point = metric_vec_orig(metric_indices);
              % generate mutants            
              for k = 1:n_seeds                    
                  % sample rates      
                  quad_exit_flag = false;
                  while ~quad_exit_flag
                      if equilibrium_flag
                          exit_flag_eq = false;
                          while ~exit_flag_eq                                
                              prop_factors = randsample(proposal_vec,8,true);                        
                              new_rate_vec = prop_factors.*rates_orig;                        
                              flux = prod(new_rate_vec(1:4))/prod(new_rate_vec(5:8));
                              new_rate_vec(5:8) = flux^.25 * new_rate_vec(5:8);
                              exit_flag_eq = all(new_rate_vec<=rate_bounds(2))&&all(new_rate_vec>=rate_bounds(1));
                          end              
                      else
                          prop_factors = randsample(proposal_vec,8,true);                        
                          new_rate_vec = prop_factors.*rates_orig;  
                          new_rate_vec(new_rate_vec>rate_bounds(2)) = rate_bounds(2);
                          new_rate_vec(new_rate_vec<rate_bounds(1)) = rate_bounds(1);                                                                       
                      end
                      % calculate eq performance metrics
                      metric_vec = calculateMetricsV2(new_rate_vec,c_range,frequency_flag);
                      new_point = metric_vec(metric_indices);
                      % check whether point moves in desired direction
                      quadrant_flag = all(ismember(new_point>orig_point,quadrant_ref));
                      % decide whther or not to take point
                      quad_exit_flag = max(quadrant_flag,rand()) > quadrant_tol;
                  end
                  % record                   
                  rate_array(k,:,b,nti) = new_rate_vec;                                                                           
                  metric_array(k,:,b,nti) = metric_vec;                                                                                                                            
              end  
          end  
      end        
      %toc         disp(['completed simulation ' sprintf('%03d',mci) ' of ' num2str(n_models)...
%             ' (model ' sprintf('%02d',mdi) ' of ' num2str(size(m_space,1)) ')'])        

  end       
  simulation_results(mdi).metric_array = metric_array;        
  simulation_results(mdi).rate_array = rate_array;    

toc
% save
fpath = [outPath save_name];
save(fpath,'simulation_results','-v7.3')