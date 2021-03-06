function [simInfo, simResults] = param_sweep_multi(metric_indices,varargin)

useParpool = 0;

% specify sampling hyperparameters
n_sim = 5; % number of independent runs to undertake
n_iters_max = 2e1; % iterations per sim
convergence_threshold = 1e-3;
nStates = 4;
testingFlag = 0 ;
numCalcFlag = 0;
sameOnRateFlag = 0;
numTesting = 0;

prop_sigma = 0.25;
n_seeds = 5; % number of new "seed" networks to generate per point on boundary
max_grid_res = 50;
min_points_per_bin = 10;
c_val = 0.95;

specFactor = 100; % ratio of specific to non-specific unbinding rates
wrongFactorConcentration = 1;

equilibrium_flag = 0; % if 1, detail balanced is enforced
half_max_flag = 0;
rnd_seed = 'shuffle';

% estimate a reasonable number of workers
myCluster = parcluster('local');
NumWorkers = ceil(myCluster.NumWorkers/2);

% names of metric options
[~,metric_names,metric_ub_vec,metric_lb_vec] = calculateMetricsMultiState([]);

%% %%%%%%%%%%%%%%%%%%%%%% Check for optional inputs %%%%%%%%%%%%%%%%%%%%%%%
for i = 1:numel(varargin)
   if ischar(varargin{i}) && i <  numel(varargin)
       eval([varargin{i} ' = varargin{i+1};'])
   end
end

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

% rate_bounds = log10([1e-4, 1e4]); % constrain transition rate magnitude
rate_bounds = repmat([-5 ; 4 ], 1, 3*nStates-4);
if size(rate_bounds,2)~=3*nStates-4
    error('inconsistent "nStates" and "rate_bounds" sizes')
end
  
% iterate through specified model specs
simInfo = struct;
simResults = struct;

% initialize random number generator
rng(rnd_seed);  

% calculate useful quantities
simInfo.rep_size = 4*max_grid_res*n_seeds;

% record relevant hyperparameters  
simInfo.c_val = c_val;
simInfo.nStates = nStates;
simInfo.rate_bounds = rate_bounds;     
simInfo.metric_names = metric_names; 
simInfo.edge_metric_indices = metric_indices;  
simInfo.metric_lb_vec = metric_lb_vec;
simInfo.metric_ub_vec = metric_ub_vec;
simInfo.equilibrium_flag = equilibrium_flag;
simInfo.specFactor = specFactor;
simInfo.wrongFactorConcentration = wrongFactorConcentration;
simInfo.half_max_flag = half_max_flag;
simInfo.testingFlag = testingFlag;
simInfo.numCalcFlag = numCalcFlag;
simInfo.sameOnRateFlag = sameOnRateFlag;
simInfo.numTesting = 0;
% set defaults and check for custom options
% simInfo = determineSweepOptions(metric_indices,varargin);
if ~simInfo.numCalcFlag
    if simInfo.half_max_flag
        % load hm functions
        loadPath = '../utilities/metricFunctions/n4_AND/';
        load([loadPath 'hmStruct.mat'],'hmStruct');
        simInfo.hmStruct = hmStruct;
    end

    simInfo = getBindingEdges(simInfo);               
end

% iterate through specified number of independent runs   
for nti = 1:n_sim
  
    %% %%%%%%%%%%%%%%%%%%%%%% Draw initial samples %%%%%%%%%%%%%%%%%%%%%%%%%%
    rate_array = NaN(n_iters_max*simInfo.rep_size+simInfo.rep_size*min_points_per_bin,3*nStates-4);
    metric_array = NaN(n_iters_max*simInfo.rep_size+simInfo.rep_size*min_points_per_bin,length(metric_names));
    area_vec = NaN(n_iters_max,1);

    % generate initial array of samples
    lb_array = repmat(rate_bounds(1,:)/4,simInfo.rep_size*min_points_per_bin,1);
    ub_array = repmat(rate_bounds(2,:)/4,simInfo.rep_size*min_points_per_bin,1);

    % call rate sampling funtion
    new_rates = sample_rates_multi(lb_array,ub_array,[],4,simInfo);
    % record
    rate_array(1:simInfo.rep_size*min_points_per_bin,:) = new_rates;

    % calculate metrics  
    if ~simInfo.numCalcFlag
        metric_array(1:simInfo.rep_size*min_points_per_bin,:) = calculateMetricsMultiState(...
                                            rate_array(1:simInfo.rep_size*min_points_per_bin,:),c_val,simInfo);
    else
        metric_array(1:simInfo.rep_size*min_points_per_bin,:) = calculateMetricsMultiStateNumSS(...
                                            rate_array(1:simInfo.rep_size*min_points_per_bin,:),c_val,simInfo);
    end
    %% %%%%%%% Perform edge-sampling until convergence or max iterations %%%%

    % initialize convergence metric
    prev_ratio = Inf;
    i_pass = 1;  
    % now perform iterative edge-biased sampling
    while prev_ratio > convergence_threshold && i_pass <= n_iters_max

        last_index = (i_pass-1)*simInfo.rep_size + simInfo.rep_size*min_points_per_bin;

        % extract and reshape arrays
        rate_array_curr = rate_array(1:last_index,:);
        metric_array_curr = metric_array(1:last_index,metric_indices);    
        use_indices = find(max(isnan(metric_array_curr),[],2)==0);

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
       
        try
            boundaryPoints = boundaryPoints(~isnan(boundaryPoints));
            sample_indices = randsample(use_indices(boundaryPoints(:)),simInfo.rep_size/n_seeds,true);
        catch
            error('afsa')
        end

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
%         n_expected = last_index+simInfo.rep_size - last_index;
%         n_actual = size(rate_array_curr(,1)*n_seeds;
        orig_rate_array = log10(repmat(rate_array_curr(sample_indices,:),n_seeds,1));

        % generate variants
        lb_array = (rate_bounds(1,:)-orig_rate_array)/prop_sigma;
        ub_array = (rate_bounds(2,:)-orig_rate_array)/prop_sigma;

        % call rate sampling funtion
        new_rates = sample_rates_multi(lb_array,ub_array,orig_rate_array,prop_sigma,simInfo);  
                                                
        % record   
        rate_array(last_index+1:last_index+simInfo.rep_size,:) = new_rates;

        % calculate metrics    
        if ~simInfo.numCalcFlag
            metric_array(last_index+1:last_index+simInfo.rep_size,:) = calculateMetricsMultiState(rate_array(last_index+1:last_index+simInfo.rep_size,:),c_val, simInfo);
        else
            metric_array(last_index+1:last_index+simInfo.rep_size,:) = calculateMetricsMultiStateNumSS(rate_array(last_index+1:last_index+simInfo.rep_size,:),c_val, simInfo);
        end
        % update metrics 
        if i_pass > 2            
            prev_ratio = area_vec(i_pass)/area_vec(i_pass-2) - 1;
        end
        i_pass = i_pass + 1;    
    end  
    % calculate cutoff point
    last_index = (i_pass-1)*simInfo.rep_size + simInfo.rep_size*min_points_per_bin;

    simResults(nti).metric_array = metric_array(1:last_index,:);
    simResults(nti).rate_array = rate_array(1:last_index,:);
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
simInfo.rate_bounds = rate_bounds;
simInfo.prop_sigma = prop_sigma;
simInfo.n_seeds = n_seeds;
simInfo.max_grid_res = max_grid_res;
simInfo.min_points_per_bin = min_points_per_bin;
simInfo.cycleTime = [];
simInfo.c_val = 1;
