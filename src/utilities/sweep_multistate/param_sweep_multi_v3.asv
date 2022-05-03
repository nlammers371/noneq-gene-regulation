function [simInfo, simResultsOut, simResultsShort] = param_sweep_multi_v3(metric_indices,functionPath,varargin)

if isempty(metric_indices)
    error('Nor metric indices specified.')
end

useParpool = 0;

% specify sampling hyperparameters
simType = 'General';
% saturationFlag = 0;
downsample_output = 0;
% testFlag = 0;
r_target = 100;
a1 = 0.98; % upper limit on induction curve region to consider
a0 = 0.02; % lower limit on induction curve region to consider
n_sim = 5; % number of independent runs to undertake
n_iters_max = 2e2; % iterations per sim
convergence_threshold = 1e-3;
ds_factor = 4e3;
PhiLimit = Inf; % entropy rate limit
TauCycleTime = 1;
% modelCollapseOption = '';
numerical_precision = 15;
numCalcFlag = 0;
% eqBindingFlag = 0;
prop_sigma = 0.25;
n_seeds = 5; % number of new "seed" networks to generate per point on boundary
max_grid_res = 50;
min_points_per_bin = 10;
cr0 = 0.95;
cr1 = 1.05;
crs = 1;
specFactor = 100; % ratio of specific to non-specific unbinding rates
cw = 1;
minError = 0.32;
% twoSpecFlag = 0;
equilibrium_flag = 0; % if 1, detail balanced is enforced
half_max_flag = 0;
rnd_seed = 'shuffle';

% estimate a reasonable number of workers
myCluster = parcluster('local');
NumWorkers = ceil(myCluster.NumWorkers);

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
  
if useParpool
    pctRunOnAll warning('off','MATLAB:rmpath:DirNotFound'); % suppress uniformative warnings
else
    warning('off','MATLAB:rmpath:DirNotFound');
end    
% pctRunOnAll warning('off','MATLAB:nearlySingularMatrix');
% iterate through specified model specs
simInfo = getSystemInfo;
simResults = struct;

% initialize random number generator
rng(rnd_seed);  

% calculate useful quantities
simInfo.rep_size = 4*max_grid_res*n_seeds;

% record relevant hyperparameters
simInfo.ds_factor = ds_factor;
simInfo.downsample_output = downsample_output;
simInfo.TauCycleTime = TauCycleTime;
% simInfo.modelCollapseOption = modelCollapseOption;
simInfo.minError = minError;
simInfo.crs = crs;
simInfo.cr1 = cr1;
simInfo.cr0 = cr0;
simInfo.a1 = a1;
simInfo.a0 = a0;
simInfo.nStates = length(simInfo.activeStateFilter);   
simInfo.r_target = r_target; % only used for Poisson noise calcualations
simInfo.edge_metric_indices = metric_indices;  
simInfo.equilibrium_flag = equilibrium_flag;
simInfo.specFactor = specFactor;
simInfo.cw = cw;
simInfo.half_max_flag = half_max_flag;
simInfo.functionPath = functionPath;
simInfo.numCalcFlag = numCalcFlag;
simInfo.numerical_precision = numerical_precision;
% simInfo.eqBindingFlag = eqBindingFlag;
simInfo.PhiLimit = PhiLimit;
% simInfo.cycleTime = cycleTime;
% simInfo.saturationFlag = saturationFlag;
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

% names of metric options
if simInfo.numCalcFlag
    [~,~,metric_names,metric_ub_vec,metric_lb_vec] = calculateMetricsNumeric_v3([]);
else    
    [~,~,metric_names,metric_ub_vec,metric_lb_vec] = calculateMetricsSym_v2([]);
end
% store
simInfo.metric_names = metric_names; 
simInfo.metric_lb_vec = metric_lb_vec;
simInfo.metric_ub_vec = metric_ub_vec;

defaultValues = simInfo.defaultValues; 
sweepVarList = simInfo.sweepVarList;

metric_names_sweep = metric_names(metric_indices);
if any(strcmp(metric_names_sweep,'CW'))
    sweepFlags(simInfo.cw_index) = true;
    try
        sweepFlags(simInfo.a_index) = false;
    catch
        sweepFlags(simInfo.b_index) = false;
    end
    simInfo.sweepFlags = sweepFlags;
    paramBounds(:,simInfo.cw_index) = [log10(specFactor/1000) ; log10(specFactor^3.1)];
    simInfo.paramBounds = paramBounds;
end

% iterate through specified number of independent runs

parfor nti = 1:n_sim%(nti = 1:n_sim,useParpool*NumWorkers) %NL: note that if useParpool=0, this will switch to serial execution for % %% %%

    %% %%%%%%%%%%%%%%%%%%%%%% Draw initial samples %%%%%%%%%%%%%%%%%%%%%%%%%%
    % initialize
    param_array = NaN(n_iters_max*simInfo.rep_size+simInfo.rep_size*min_points_per_bin,length(sweepVarList));
    metric_array = NaN(n_iters_max*simInfo.rep_size+simInfo.rep_size*min_points_per_bin,length(metric_names));
    iter_id_vec = NaN(n_iters_max*simInfo.rep_size+simInfo.rep_size*min_points_per_bin,1);
    area_vec = NaN(n_iters_max,1);

    % set param values for static variables
    param_array(:,~sweepFlags) = repmat(defaultValues(~sweepFlags),size(param_array,1),1);
    
    % generate initial array of samples
    lb_array = repmat(paramBounds(1,sweepFlags)/4,simInfo.rep_size*min_points_per_bin,1);
    ub_array = repmat(paramBounds(2,sweepFlags)/4,simInfo.rep_size*min_points_per_bin,1);

    % call rate sampling funtion
    temp_slice = param_array(1:simInfo.rep_size*min_points_per_bin,:);
    temp_slice(:,sweepFlags) = sample_rates_multi_v2(lb_array,ub_array,[],4,simInfo);        
    iter_id_vec(1:simInfo.rep_size*min_points_per_bin) = 1;
    
    % calculate metrics  
    if simInfo.numCalcFlag
        [metric_array(1:simInfo.rep_size*min_points_per_bin,:), temp_slice_adjusted]...
                                                = calculateMetricsNumeric_v3(...
                                                  temp_slice,simInfo);
    else
        [metric_array(1:simInfo.rep_size*min_points_per_bin,:), temp_slice_adjusted]...
                                                = calculateMetricsSym_v2(...
                                                  temp_slice,simInfo);
    end
    
    % record
    param_array(1:simInfo.rep_size*min_points_per_bin,:) = temp_slice_adjusted;
    
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
          
            [boundaryPoints, metric_array_filtered] = findBoundaryPoints(...
                          metric_array_curr,use_indices,min_points_per_bin,max_grid_res);

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

        temp_slice = param_array(last_index+1:last_index+simInfo.rep_size,:);
        iter_id_vec(last_index+1:last_index+simInfo.rep_size) = i_pass+1;
        temp_slice(:,sweepFlags) = new_params;
        % calculate metrics    
        if simInfo.numCalcFlag
            [metric_array(last_index+1:last_index+simInfo.rep_size,:),temp_slice_adjusted]...
                        = calculateMetricsNumeric_v3(temp_slice, simInfo);
        else
            [metric_array(last_index+1:last_index+simInfo.rep_size,:),temp_slice_adjusted]...
                        = calculateMetricsSym_v2(temp_slice, simInfo);
        end
        
        % record   
        param_array(last_index+1:last_index+simInfo.rep_size,:) = temp_slice_adjusted;
        
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
    simResults(nti).iter_id_vec = iter_id_vec(1:last_index,:);
    simResults(nti).sim_id_vec = nti*ones(size(iter_id_vec));
    simResults(nti).area_vec = area_vec(1:i_pass-1);
    simResults(nti).n_iters = i_pass-1;
    simResults(nti).convergence_flag = i_pass-1 < n_iters_max;
end       

%% %%%%%%%%%%%%%%%%%%%%%% Update info structure %%%%%%%%%%%%%%%%%%%%%%%%%%%
   
simInfo.half_max_flag = half_max_flag;
simInfo.equilibrium_flag = equilibrium_flag;
simInfo.rnd_seed = rnd_seed;
simInfo.n_sim = n_sim;
simInfo.n_iters_max = n_iters_max;
simInfo.paramBounds = paramBounds;
simInfo.prop_sigma = prop_sigma;
simInfo.n_seeds = n_seeds;
simInfo.max_grid_res = max_grid_res;
simInfo.min_points_per_bin = min_points_per_bin;
simInfo.c_val = crs;

% add diagnostic fields from sim struct
if n_sim > 0
    simInfo.convergence_flag_vec = [simResults.convergence_flag];
    simInfo.n_iter_vec = [simResults.n_iters];

    % concatenate key sim struct fields
    simResultsOut = struct;
    simResultsOut.metric_array = vertcat(simResults.metric_array);
    simResultsOut.rate_array = vertcat(simResults.rate_array);
    simResultsOut.iter_id_vec = vertcat(simResults.iter_id_vec);
    simResultsOut.sim_id_vec = vertcat(simResults.sim_id_vec);
    simResultsOut.areaArray = vertcat(simResults.area_vec);

    % remove NaN values
    nan_ft = any(isnan(simResultsOut.rate_array),2) | any(isnan(simResultsOut.metric_array(:,metric_indices)),2);
    simResultsOut.metric_array = simResultsOut.metric_array(~nan_ft,:);
    simResultsOut.rate_array = simResultsOut.rate_array(~nan_ft,:);
    simResultsOut.iter_id_vec = simResultsOut.iter_id_vec(~nan_ft);
    simResultsOut.sim_id_vec = simResultsOut.sim_id_vec(~nan_ft);
    
    if simInfo.downsample_output
        % identify and keep global edge points
        use_indices = 1:size(simResultsOut.metric_array,1);    
        [boundaryPoints, ~,randomPoints] = findBoundaryPoints(...
                              simResultsOut.metric_array(:,metric_indices),...
                              use_indices,min_points_per_bin,3e3,ds_factor,true);

        % randomly select 1% of other datapoints to keep
        simResultsShort = simResultsOut;
        keepPoints = [boundaryPoints ; randomPoints];
        simResultsShort.metric_array = simResultsShort.metric_array(keepPoints,:);
        simResultsShort.rate_array = simResultsShort.rate_array(keepPoints,:);
    %     simResultsShort
    else
        simResultsShort = [];
    end   
else    
    simResultsShort = [];
    simResultsOut = [];
end    