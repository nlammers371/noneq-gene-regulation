clear 
close all

addpath('../utilities')
outPath = '../../out/param_sweep/';
mkdir(outPath)

%%%% Set parameter defaults

% define model characteristics
m_space = [1 7; 3 5];
activator_flag_vec = [1 0];
edge_key = [7 8 5 6 3 4 1 2];
metric_indices = 1:3;
metric_names = {'fidelity','dynamicRange','info_rate'};
opt_options = optimoptions('fmincon','Display','off');
% generate all possible group sizes
base = 1:4;
edge_samp_cell = {base};
g_size = 1:4;
edge_options = 5:8;
for g = g_size
    c_array = nchoosek(edge_options,g);
    for c_vec = 1:size(c_array,1)
        edge_samp_cell = [edge_samp_cell{:} {[base c_array(c_vec,:)]}];
    end
end

% specify sampling hyperparameters
n_reps =  5; % number of initializations for optimization routine
c_vec = [.9 1.1];
rate_bounds = [1e-4, 1e2]; % constrain transition rate magnitude


% calculate useful quantities
n_models = numel(edge_samp_cell);
edge_vec = 1:8;
% names of metric options

% set save name
save_name = 'param_point_optimization.mat';

% iterate through specified model specs
opt_results = struct;
for mdi = 1:size(m_space,1)  
    activator_flag = activator_flag_vec(mdi);
    % initialize arrays
    rate_array = NaN(n_reps,8,numel(metric_names),n_models);  
    metric_array = NaN(n_reps,numel(metric_names),n_models);               
    % record relevant hyperparameters        
    opt_results(mdi).model_spec = m_space(mdi,:);                    
    opt_results(mdi).c_vec = c_vec; 
    opt_results(mdi).rate_bounds = rate_bounds; 
    opt_results(mdi).edge_samp_cell = edge_samp_cell; 
    opt_results(mdi).metric_names = metric_names; 
    opt_results(mdi).edge_metric_indices = metric_indices;     
    % iterate through specified number of independent MCMC runs   
    for mci = 1:n_models    
        tic
        opt_edges = edge_samp_cell{mci};
        n_samp_edges = numel(opt_edges);
        sym_edges = edge_vec(~ismember(edge_vec,opt_edges));              
        for m = 1:numel(metric_indices)
            tic
            metric_index = metric_indices(m);
            call_opt_fun = @(opt_rates) -rate_opt_general(opt_rates,opt_edges,metric_index,c_vec,activator_flag);
            parfor b = 1:n_reps
                % decide on initial "move"            
                init_rate_vec = NaN(1,8);
                init_rates = randsample(logspace(log10(rate_bounds(1)),log10(rate_bounds(2))),numel(opt_edges),true);                                      
                % perform optimization
                rates_fit = fmincon(call_opt_fun,init_rates,[],[],[],[],repelem(rate_bounds(1),numel(opt_edges)),...
                    repelem(rate_bounds(2),numel(opt_edges)),[],opt_options);
                rate_vec_out = NaN(1,8);
                rate_vec_out(opt_edges) = rates_fit;
                rate_vec_out(sym_edges) = rate_vec_out(edge_key(sym_edges));                
                rate_array(b,:,m,mci) = rate_vec_out;                 
                % calculate  metrics
                metric_array(b,m,mci) = -call_opt_fun(rates_fit);                                                    
            end
            toc
        end
    end
    opt_results(mdi).metric_array = metric_array;        
    opt_results(mdi).rate_array = rate_array;    
end
% save
fpath = [outPath save_name];
save(fpath,'opt_results','-v7.3')