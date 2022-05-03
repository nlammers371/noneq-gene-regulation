function [simulation_results, fpath] = param_sweep(metric_indices,varargin)
addpath('../utilities')
outPath = '../../out/param_sweep/';
mkdir(outPath)

%%%% Set parameter defaults

% define model characteristics
m_space = [1 7; 3 5];
activator_flag_vec = [1 0];
edge_key = [7 8 5 6 3 4 1 2];
edge_samp_cell = {1:4, 1:6, 1:7, [1:6 8], 1:8};
% edge_samp_cell = {1:8};

% specify sampling hyperparameters
n_sim = 5; % number of independent runs to undertake
n_iters = 2e1; % iterations per sim
c_range = linspace(2/11,20/11,101);
rate_bounds = [1e-4, 1e1]; % constrain transition rate magnitude
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
edge_vec = 1:8;
n_models = numel(edge_samp_cell);
n_samples_init = n_seeds*n_bound_samples;
mem_steps = ceil(n_b_points_max/n_samples_init);
% names of metric options
metric_names = {'fold-concentration','fold-production',...
    'fidelity','dynamicRange','dmdc','precision','decision_time','error_rate'};
% set save name
save_name = ['param_sweep_results_' metric_names{metric_indices(1)} '_' ...
    metric_names{metric_indices(2)} '.mat'];

% iterate through specified model specs
simulation_results = struct;
for mdi = 1:size(m_space,1)  
    activator_flag = activator_flag_vec(mdi);
    % initialize arrays
    rate_array = NaN(n_seeds,8,n_bound_samples,n_iters*n_sim,n_models);  
    metric_array = NaN(n_seeds,numel(metric_names),n_bound_samples,n_iters*n_sim,n_models);               
    % record relevant hyperparameters        
    simulation_results(mdi).model_spec = m_space(mdi,:);                    
    simulation_results(mdi).c_range = c_range; 
    simulation_results(mdi).rate_bounds = rate_bounds; 
    simulation_results(mdi).edge_samp_cell = edge_samp_cell; 
    simulation_results(mdi).metric_names = metric_names; 
    simulation_results(mdi).edge_metric_indices = metric_indices;     
    % iterate through specified number of independent MCMC runs   
    for mci = 1:n_models    
        tic
        sample_edges = edge_samp_cell{mci};
        n_samp_edges = numel(sample_edges);
        sym_edges = edge_vec(~ismember(edge_vec,sample_edges));
        restart_indices = 1:n_iters:n_sim*n_iters;
        for nti = 1:n_iters*n_sim
            % if at a restart index, draw set of initial points
            if ismember(nti,restart_indices)
                for b = 1:n_bound_samples
                    for s = 1:n_seeds
                        % decide on initial "move"            
                        orig_rate_vec = NaN(1,8);
                        orig_rate_vec(1:4) = randsample(logspace(log10(rate_bounds(1)),log10(rate_bounds(2))),4,true);                       
                        if all(ismember(sample_edges,1:4))
                            orig_rate_vec(sym_edges) = orig_rate_vec(edge_key(sym_edges));
                        else
                            % check which edges lie opposite of constrained
                            % edges
                            opp_edges = sym_edges-2;                            
                            % assign free edges
                            free_edges = edge_vec(~ismember(edge_vec,[1:4, sym_edges, opp_edges]));
                            % sample opposite edges subject to adjusted
                            % bounds
                            for i = 1:numel(opp_edges)
                                p_rate = orig_rate_vec(opp_edges(i)-4);
                                e_rate = orig_rate_vec(opp_edges(i)-2);
                                ub = min((e_rate*p_rate)/rate_bounds(1),rate_bounds(2));
                                lb = max((e_rate*p_rate)/rate_bounds(2),rate_bounds(1));
                                orig_rate_vec(opp_edges(i)) = randsample(logspace(log10(lb),log10(ub)),1,true);                                              
                            end                            
                            orig_rate_vec(free_edges) = randsample(logspace(log10(rate_bounds(1)),log10(rate_bounds(2))),numel(free_edges),true);                       
                            % assign determined edges
                            orig_rate_vec(sym_edges) = orig_rate_vec(sym_edges-4) .* orig_rate_vec(sym_edges-6) ./ orig_rate_vec(sym_edges-2);
                        end
                        rate_array(s,:,b,nti,mci) = orig_rate_vec;                                                                
                        % calculate  metrics
                        metric_array(s,:,b,nti,mci) = calculateMetricsConstrained(orig_rate_vec,c_range,activator_flag);                                        
                    end
                end  
            else                 
                % calculate lookback distance
                mem_max = mod(nti-1,n_iters);
                mem = min(mem_max,mem_steps);
                % extract and reshape arrays
                rate_array_curr = reshape(permute(rate_array(:,:,:,nti-mem:nti-1,mci),[1 4 3 2]),[],8);
                metric_array_curr = reshape(permute(metric_array(:,metric_indices,:,nti-mem:nti-1,mci),[1 4 3 2]),[],2);    
                use_indices = 1:size(metric_array_curr,1);
                if any(metric_indices==5)
                    use_indices = find(metric_array_curr(:,metric_indices==5)>.01);
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
                    % generate mutants            
                    for k = 1:n_seeds                    
                        % sample rates
                        prop_factors = randsample(proposal_vec,n_samp_edges,true);
                        new_rate_vec = NaN(1,8);
                        new_rate_vec(sample_edges) = prop_factors.*rates_orig(sample_edges);                        
                        new_rate_vec(new_rate_vec>rate_bounds(2)) = rate_bounds(2);
                        new_rate_vec(new_rate_vec<rate_bounds(1)) = rate_bounds(1);
                        if all(ismember(sample_edges,1:4))
                            new_rate_vec(sym_edges) = new_rate_vec(edge_key(sym_edges));
                        else
                            % check which edges lie opposite of constrained
                            % edges
                            opp_edges = sym_edges-2;  
                            % impose limits
                            for i = 1:numel(opp_edges)
                                p_rate = new_rate_vec(opp_edges(i)-4);
                                e_rate = new_rate_vec(opp_edges(i)-2);
                                % impose bounds
                                new_rate_vec(opp_edges(i)) = min((e_rate*p_rate)/rate_bounds(1),...
                                    new_rate_vec(opp_edges(i)));
                                new_rate_vec(opp_edges(i)) = max((e_rate*p_rate)/rate_bounds(2),...
                                    new_rate_vec(opp_edges(i)));                                
                            end
                            % assign free edges
%                             free_edges = edge_vec(~ismember(edge_vec,[1:4, sym_edges, opp_edges]));
%                             new_rate_vec(free_edges) = randsample(logspace(log10(rate_bounds(1)),log10(rate_bounds(2))),numel(free_edges),true);                       
                            % assign determined edges
                            new_rate_vec(sym_edges) = new_rate_vec(sym_edges-4) .* new_rate_vec(sym_edges-6) ./ new_rate_vec(sym_edges-2);
                        end                                                
                        rate_array(k,:,b,nti,mci) = new_rate_vec;               
                        % calculate eq performance metrics
                        metric_vec = calculateMetricsConstrained(new_rate_vec,c_range,activator_flag);                                   
                        metric_array(k,:,b,nti,mci) = metric_vec;                                                                                                                            
                    end  
                end  
            end
        end
        disp(['completed simulation ' sprintf('%03d',mci) ' of ' num2str(n_models)...
            ' (model ' sprintf('%02d',mdi) ' of ' num2str(size(m_space,1)) ')'])        
        toc
    end       
    simulation_results(mdi).metric_array = metric_array;        
    simulation_results(mdi).rate_array = rate_array;    
end
% save
fpath = [outPath save_name];
save(fpath,'simulation_results','-v7.3')