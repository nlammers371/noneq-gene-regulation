function [simInfo, simResults] = determineSweepOptions(metric_indices,varargin)

%% %%%%%%%%%%%%%%%%%%%%%% Set parameter defaults %%%%%%%%%%%%%%%%%%%%%%%%%%

useParpool = 0;

% specify sampling hyperparameters
n_sim = 5; % number of independent runs to undertake
% n_iters_max = 2e1; % iterations per sim
% convergence_threshold = 1e-3;
nStates = 4;
testingFlag = 0 ;
numTesting = 0;
sameOnRateFlag = 0;

% prop_sigma = 0.25;
n_seeds = 5; % number of new "seed" networks to generate per point on boundary
max_grid_res = 50;
% min_points_per_bin = 10;
c_val = 0.95;

specFactor = 100; % ratio of specific to non-specific unbinding rates
wrongFactorConcentration = 1;

equilibrium_flag = 0; % if 1, detail balanced is enforced
half_max_flag = 0;
rnd_seed = 'shuffle';

% estimate a reasonable number of workers
myCluster = parcluster('local');
NumWorkers = ceil(myCluster.NumWorkers/3);

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
% calculate useful quantities
rep_size = 4*max_grid_res*n_seeds;
  
% iterate through specified model specs
simInfo = struct;
simResults = struct;
% initialize random number generator
rng(rnd_seed);  
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
simInfo.numTesting = numTesting;
simInfo.sameOnRateFlag = sameOnRateFlag;