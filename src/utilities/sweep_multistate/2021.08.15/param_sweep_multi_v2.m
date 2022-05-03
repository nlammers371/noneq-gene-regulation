function [simInfo, simResults] = param_sweep_multi_v2(metric_indices,functionPath,varargin)

useParpool = 0;

% specify sampling hyperparameters
simType = 'General';
saturationFlag = 0;
testFlag = 0;
n_sim = 5; % number of independent runs to undertake
n_iters_max = 2e1; % iterations per sim
convergence_threshold = 1e-3;
PhiLimit = Inf; % entropy rate limit
TauCycleLimit = 0;
nStates = 4;
numerical_precision = 15;
numCalcFlag = 0;
eqBindingFlag = 0;
prop_sigma = 0.25;
n_seeds = 5; % number of new "seed" networks to generate per point on boundary
max_grid_res = 25;
min_points_per_bin = 10;
cr0 = 0.95;
cr1 = 1.05;
crs = 1;
specFactor = 100; % ratio of specific to non-specific unbinding rates
cw = 1;
minError = 0.32;
twoSpecFlag = 0;
equilibrium_flag = 0; % if 1, detail balanced is enforced
half_max_flag = 0;
rnd_seed = 'shuffle';

% estimate a reasonable number of workers
myCluster = parcluster('local');
NumWorkers = ceil(myCluster.NumWorkers/2);

%% %%%%%%%%%%%%%%%%%%%%%% Check for optional inputs %%%%%%%%%%%%%%%%%%%%%%%
for i = 1:numel(varargin)
   if ischar(varargin{i}) && i <  numel(varargin)
       eval([varargin{i} ' = varargin{i+1};'])
   end
end

if ~exist(functionPath,'dir')
    error('Invalid metric function path provided. Check second input argument')
end

% make sure we're linked to the appropriate function subfolder
warning('off','MATLAB:rmpath:DirNotFound') % suppress uniformative warnings
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));

%% %%%%%%%%%%%%%%%%% calculate params and initialize arrays %%%%%%%%%%%%%%%
NumWorkers = min([n_sim, NumWorkers]);

% initialize parpool
if useParpool
    p = gcp('nocreate');
    if isempty(p)
      parpool(NumWorkers);
    elseif p.NumWorkers~=NumWorkers
      delete(p);
      parpool(NumWorkers);
    end
end
  
if nStates == 18 % no symbolic option for 18 states...
    numCalcFlag = 1;
end

% names of metric options
if numCalcFlag
    [~,metric_names,metric_ub_vec,metric_lb_vec] = calculateMetricsNumeric_v3([]);
else    
    [~,metric_names,metric_ub_vec,metric_lb_vec] = calculateMetricsSym([]);
end

% iterate through specified model specs
simInfo = getSystemInfo;
simResults = struct;

% initialize random number generator
rng(rnd_seed);  

% calculate useful quantities
simInfo.rep_size = 4*max_grid_res*n_seeds;

% record relevant hyperparameters
simInfo.minError = minError;
simInfo.crs = crs;
simInfo.cr1 = cr1;
simInfo.cr0 = cr0;
simInfo.nStates = nStates;   
simInfo.metric_names = metric_names; 
simInfo.edge_metric_indices = metric_indices;  
simInfo.metric_lb_vec = metric_lb_vec;
simInfo.metric_ub_vec = metric_ub_vec;
simInfo.equilibrium_flag = equilibrium_flag;
simInfo.specFactor = specFactor;
simInfo.cw = cw;
simInfo.half_max_flag = half_max_flag;
simInfo.functionPath = functionPath;
simInfo.numCalcFlag = numCalcFlag;
simInfo.numerical_precision = numerical_precision;
simInfo.eqBindingFlag = eqBindingFlag;
simInfo.PhiLimit = PhiLimit;
simInfo.TauCycleLimit = TauCycleLimit;
simInfo.testFlag = testFlag;
simInfo.saturationFlag = saturationFlag;
simInfo.convergence_threshold = convergence_threshold;

% override default bounds if option was passed 
if exist('paramBounds','var')
    simInfo.paramBounds = paramBounds;
end
if exist('sweepFlags','var')
    simInfo.sweepFlags = sweepFlags;
end
      
% deal with Half-max options
simInfo = processHMOptions(simInfo);

% extract a few useful variables into workspace
paramBounds = simInfo.paramBounds;
sweepFlags = simInfo.sweepFlags;

% set parameter bound and value options based on simType
simInfo = set_parameter_options(simInfo,simType);

defaultValues = simInfo.defaultValues; 
sweepVarList = simInfo.sweepVarList;


metric_names_sweep = metric_names(metric_indices);
if any(strcmp(metric_names_sweep,'CW'))
    sweepFlags(simInfo.cw_index) = true;
    simInfo.sweepFlags = sweepFlags;
    paramBounds(:,simInfo.cw_index) = [0 ; log10(specFactor^3.1)];
    simInfo.paramBounds = paramBounds;
end

% iterate through specified number of independent runs

for nti = 1:n_sim%parfor (nti = 1:n_sim,useParpool*NumWorkers) %NL: note that if useParpool=0, this will switch to serial execution

    %% %%%%%%%%%%%%%%%%%%%%%% Draw initial samples %%%%%%%%%%%%%%%%%%%%%%%%%%
    % initialize
    param_array = NaN(n_iters_max*simInfo.rep_size+simInfo.rep_size*min_points_per_bin,length(sweepVarList));
    metric_array = NaN(n_iters_max*simInfo.rep_size+simInfo.rep_size*min_points_per_bin,length(metric_names));
    area_vec = NaN(n_iters_max,1);

    % set param values for static variables
    param_array(:,~sweepFlags) = repmat(defaultValues(~sweepFlags),size(param_array,1),1);
    
    % generate initial array of samples
    lb_array = repmat(paramBounds(1,sweepFlags)/4,simInfo.rep_size*min_points_per_bin,1);
    ub_array = repmat(paramBounds(2,sweepFlags)/4,simInfo.rep_size*min_points_per_bin,1);

    % call rate sampling funtion
    new_rates = sample_rates_multi_v2(lb_array,ub_array,[],4,simInfo);
    
    % record
    param_array(1:simInfo.rep_size*min_points_per_bin,sweepFlags) = new_rates;

    % calculate metrics  
    if simInfo.numCalcFlag
        [metric_array(1:simInfo.rep_size*min_points_per_bin,:)]...
                                                = calculateMetricsNumeric_v3(...
                                                param_array(1:simInfo.rep_size*min_points_per_bin,:),simInfo);
    else
        [metric_array(1:simInfo.rep_size*min_points_per_bin,:)]...
                                                = calculateMetricsSym(...
                                                param_array(1:simInfo.rep_size*min_points_per_bin,:),simInfo);
    end
    %% %%%%%%% Perform edge-sampling until convergence or max iterations %%%%

    % initialize convergence metric
    prev_ratio = Inf;
    i_pass = 1;  
    % now perform iterative edge-biased sampling
    while prev_ratio > convergence_threshold && i_pass <= n_iters_max

        last_index = (i_pass-1)*simInfo.rep_size + simInfo.rep_size*min_points_per_bin;

        % extract and reshape arrays
        rate_array_curr = param_array(1:last_index,:);
        metric_array_curr = metric_array(1:last_index,metric_indices);    
        use_indices = find(max(isnan(metric_array_curr),[],2)==0);

        if length(use_indices) > 100
            % define mesh for edge sampling
            metric_array_filtered = metric_array_curr(use_indices,:);
            nPoints = size(metric_array_filtered,1);
            nEdges = min(ceil(nPoints/min_points_per_bin),max_grid_res)+1;
            grid1 = quantile(unique(metric_array_filtered(:,1)),linspace(0,1,nEdges));
            grid2 = quantile(unique(metric_array_filtered(:,2)),linspace(0,1,nEdges));

            % assign obs to bins
            [~,~,~,bin1,bin2] = histcounts2(metric_array_filtered(:,1),metric_array_filtered(:,2),grid1,grid2);   

            % find edge values      
            EdgeIndices1 = NaN(1,2*(nEdges-1));
            for b = 1:2:2*length(grid1)
                binIndices = find(bin2==ceil(b/2));
                if ~isempty(binIndices)
                    [~, min_i] = min(metric_array_filtered(binIndices,1));
                    EdgeIndices1(b) = binIndices(min_i);  
                    [~, max_i] = max(metric_array_filtered(binIndices,1));
                    EdgeIndices1(b+1) = binIndices(max_i);  
                end
            end

            EdgeIndices2 = NaN(1,2*(nEdges-1));
            for b = 1:2:2*length(grid1)
                binIndices = find(bin1==ceil(b/2));
                if ~isempty(binIndices)
                    [~, min_i] = min(metric_array_filtered(binIndices,2));
                    EdgeIndices2(b) = binIndices(min_i);  
                    [~, max_i] = max(metric_array_filtered(binIndices,2));
                    EdgeIndices2(b+1) = binIndices(max_i);  
                end
            end

            boundaryPoints = [EdgeIndices1' EdgeIndices2'];  

            boundaryPoints = boundaryPoints(~isnan(boundaryPoints));
            sample_indices = randsample(use_indices(boundaryPoints(:)),simInfo.rep_size/n_seeds,true);     

            % calculate area enclosed by points
            xPoints = metric_array_filtered(boundaryPoints,1);
            yPoints = metric_array_filtered(boundaryPoints,2);
            [theta,~] = cart2pol(xPoints-mean(xPoints),yPoints-mean(yPoints));
            [~, si] = sort(theta);%sgolayfilt(metric_array_filtered(boundaryPoints,1), 2, 3);
            xPoints = xPoints(si);
            yPoints = yPoints(si);%sgolayfilt(metric_array_filtered(boundaryPoints,2), 2, 3);
            area_vec(i_pass) = polyarea(xPoints,yPoints);

            % generate expanded array for sampling. Dynamically increase expansion factor
            % if too many rejected samples from previous iteration
            orig_param_array = log10(repmat(rate_array_curr(sample_indices,:),n_seeds,1));

            % generate variants
            lb_array = (paramBounds(1,sweepFlags)-orig_param_array(:,sweepFlags))/prop_sigma;
            ub_array = (paramBounds(2,sweepFlags)-orig_param_array(:,sweepFlags))/prop_sigma;

            % call rate sampling funtion
            new_params = sample_rates_multi_v2(lb_array,ub_array,orig_param_array(:,sweepFlags),prop_sigma,simInfo);  
            
            % indicate that sampling worked
            err_flag = 0;
        else
            lb_array = repmat(paramBounds(1,sweepFlags)/4,simInfo.rep_size,1);
            ub_array = repmat(paramBounds(2,sweepFlags)/4,simInfo.rep_size,1);
            
            new_params = sample_rates_multi_v2(lb_array,ub_array,[],4,simInfo);
            
            area_vec(i_pass) = realmin;
            
            % indicate that sampling had to be reset
            err_flag = 1;
        end
        
        % record   
        param_array(last_index+1:last_index+simInfo.rep_size,sweepFlags) = new_params;

        % calculate metrics    
        if simInfo.numCalcFlag
            [metric_array(last_index+1:last_index+simInfo.rep_size,:)]...
                        = calculateMetricsNumeric_v3(...
                        param_array(last_index+1:last_index+simInfo.rep_size,:), simInfo);
        else
            [metric_array(last_index+1:last_index+simInfo.rep_size,:)]...
                        = calculateMetricsSym(...
                        param_array(last_index+1:last_index+simInfo.rep_size,:), simInfo);
        end
        % update metrics 
        if i_pass > 3      
            if ~err_flag
                prev_ratio1 = area_vec(i_pass)/area_vec(i_pass-2) - 1;
                prev_ratio2 = area_vec(i_pass-1)/area_vec(i_pass-3) - 1;
                prev_ratio = max([prev_ratio1,prev_ratio2]);
            else
                prev_ratio = Inf;
            end
        end
        i_pass = i_pass + 1;    
    end  
    % calculate cutoff point
    last_index = (i_pass-1)*simInfo.rep_size + simInfo.rep_size*min_points_per_bin;

    simResults(nti).metric_array = metric_array(1:last_index,:);
    simResults(nti).rate_array = param_array(1:last_index,:);
    simResults(nti).area_vec = area_vec(1:i_pass-1);
    simResults(nti).n_iters = i_pass-1;
    simResults(nti).convergence_flag = i_pass-1<n_iters_max;
end       

%% %%%%%%%%%%%%%%%%%%%%%% Update info structure %%%%%%%%%%%%%%%%%%%%%%%%%%%
   
simInfo.half_max_flag = half_max_flag;
simInfo.equilibrium_flag = equilibrium_flag;
simInfo.rnd_seed = rnd_seed;
simInfo.n_sim = n_sim;
simInfo.n_iters_max = n_iters_max;
simInfo.param_bounds = paramBounds;
simInfo.prop_sigma = prop_sigma;
simInfo.n_seeds = n_seeds;
simInfo.max_grid_res = max_grid_res;
simInfo.min_points_per_bin = min_points_per_bin;
simInfo.cycleTime = [];
simInfo.c_val = 1;
