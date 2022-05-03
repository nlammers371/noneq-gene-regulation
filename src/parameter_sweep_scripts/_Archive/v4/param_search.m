clear 
close all

addpath('../utilities')
outPath = '../../out/bivariate_parameter_sweeps/';
mkdir(outPath)

%%%% Set parameter defaults
equilibrium_flag = 0;
metric_indices = [4,6];

% define range over which to explore independent variable
var1_range = linspace(0,.5,26);
n_mcmc = 100;
% define model characteristics
m_space = [1 7; 3 5];
frequency_flag_vec = [1 0];
% specify sampling hyperparameters
c_range = linspace(2/2.1,2/2.1*1.1,100);
rate_bounds = [1e-4, 1e4]; % constrain transition rate magnitude
proposal_vec = logspace(-.5,.5,1000);


% check for optional inputs
% for i = 1:numel(varargin)
%    if ischar(varargin{i}) && i <  numel(varargin)
%        eval([varargin{i} ' = varargin{i+1};'])
%    end
% end

% names of metric options
[~, metric_names] = calculateMetricsV2([],[],[]);
% set save name
save_name = ['param_search_results_' metric_names{metric_indices(1)} '_' ...
    metric_names{metric_indices(2)} 'eq' num2str(equilibrium_flag) '.mat'];
tic
% iterate through specified model specs
simulation_results = struct;
for mdi = 1:size(m_space,1)  
    frequency_flag = frequency_flag_vec(mdi);
    % initialize arrays
    rate_array = NaN(n_mcmc,8);  
    metric_array = NaN(n_mcmc,numel(metric_names));               
    % record relevant hyperparameters        
    simulation_results(mdi).model_spec = m_space(mdi,:);                    
    simulation_results(mdi).c_range = c_range; 
    simulation_results(mdi).rate_bounds = rate_bounds;     
    simulation_results(mdi).metric_names = metric_names; 
    simulation_results(mdi).edge_metric_indices = metric_indices;     
    simulation_results(mdi).equilibrium_flag = equilibrium_flag; 
    
    % iterate through var range
    for i = 1:numel(var1_range)-1
        var1_lb = var1_range(1);
        var1_ub = var1_range(2);                
        %%%%%%%%%%%%%%%%%%%%%%%
        % decide on initial move        
        exit_flag_bound = false;
        while ~exit_flag_bound
            % impose detailed balance if necessary
            if equilibrium_flag
                exit_flag_eq = false;
                while ~exit_flag_eq
                    orig_rate_vec = randsample(logspace(log10(rate_bounds(1)),log10(rate_bounds(2))),8,true);                                              
                    flux = prod(orig_rate_vec(1:4))/prod(orig_rate_vec(5:8));
                    orig_rate_vec(5:8) = flux^.25 * orig_rate_vec(5:8);
                    exit_flag_eq = all(orig_rate_vec<=rate_bounds(2))&&all(orig_rate_vec>=rate_bounds(1));
                end                            
            else
                orig_rate_vec = randsample(logspace(log10(rate_bounds(1)),log10(rate_bounds(2))),8,true);                                              
            end
            % calculate metrics
            metric_vec = calculateMetricsV2(orig_rate_vec,c_range,frequency_flag);   
            % check that we are in correct range for param 1
            val1 = metric_vec(metric_indices(1));
            exit_flag_bound = val1>=var1_lb & val1<var1_ub;
        end                                                                      
        % record
        rate_array(1,:) = orig_rate_vec;  
        metric_array(1,:) = metric_vec;
    end
    % now attempt to optimize wrpt var2
    for m = 2:n_mcmc
        rates_orig = rate_array(m-1,:);
        metric_orig = metric_array(m-1,metric_indices(2));
        exit_flag_bound = false;
        while ~exit_flag_bound
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
            val1 = metric_vec(metric_indices(1));
            % check whether point is within desired bounds
            exit_flag_bound = val1>=var1_lb & val1<var1_ub;
        end
        % now perform MH move
        metric_new = metric_vec(metric_indices(2));
        MH_factor = min(metric_new/metric_orig,1) > rand();
        if MH_factor
            rate_array(m,:) = new_rate_vec;
            metric_array(m,:) = metric_vec;
        else
            rate_array(m,:) = rate_array(m-1,:);
            metric_array(m,:) = metric_array(m-1,:);
        end                                                                                                                                 
    end  
end