% script to look at tradeoff between H0 and f

clear 
close all
addpath(genpath('../utilities/'))

% set basic parameters
nStates = 6;
paramBounds = repmat([-10 ; 6],1,3*nStates-4); % constrain transition rate magnitude
[~,~,metric_names] = calculateMetricsSym_v2([]);

% specify function path
functionPath = ['../utilities/metricFunctions/n' num2str(nStates) '_OR/'];

% make sure we're linked to the appropriate function subfolder% make sure we're linked to the appropriate function subfolder
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));

% define save path

DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\manuscript\';
% DropboxFolder = 'S:\Nick\Dropbox\Nonequilibrium\Nick\manuscript\';

FigPath = [DropboxFolder 'performance_metrics_vs_cw' filesep];
mkdir(FigPath);         

% get index of useful metrics
spec_index = find(strcmp(metric_names,'Specificity'));
sharp_right_norm_index = find(strcmp(metric_names,'SharpnessRightNorm'));
ir_index = find(strcmp(metric_names,'DecisionRateNorm'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
precision_index = find(strcmp(metric_names,'Precision'));
rate_index = find(strcmp(metric_names,'Production Rate'));
cw_index = find(strcmp(metric_names,'CW'));

% set sim options
sweep_options = {'n_seeds',5,'n_iters_max',50,'n_sim',10,'nStates',nStates};

% calculate sensitivity bound
alpha_factor = 100;
% seq = 1/4;

% specify plot options 
markerAlpha = 0.5; % marker transparency
markerSize = 75; % marker size
n_plot = 3e3;
bit_factor = log2(exp(1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run parameter sweeps

cw = 1; % precise value unimportant...just something small enough to be negligible
tic
[sim_info_sharp, sim_struct_sharp] = param_sweep_multi_v3([cw_index ir_index],functionPath,sweep_options{:},...
                                          'half_max_flag',true,...'cw',1e1,...
                                          'equilibrium_flag',false,'specFactor',alpha_factor);
                                        
% [sim_info_spec, sim_struct_spec] = param_sweep_multi_v3([spec_index ir_index],functionPath,sweep_options{:},...
%                                           'half_max_flag',true,'cw',1e2,...
%                                           'equilibrium_flag',false,'specFactor',alpha_factor);                                        

toc
%%

% extract vectors
metric_array_sharp = sim_struct_sharp.metric_array;   
rate_array_sharp = sim_struct_sharp.rate_array;   
sharpness_vec_sharp = metric_array_sharp(:,sharpness_index);
ir_vec_sharp = metric_array_sharp(:,ir_index);
precision_vec_sharp = sqrt(exp(metric_array_sharp(:,precision_index)));
spec_vec_sharp = metric_array_sharp(:,spec_index);
pd_vec_sharp = metric_array_sharp(:,rate_index);
cw_vec_sharp = metric_array_sharp(:,cw_index);

% generate axis
cw_axis = linspace(0,4,26);

n_keep = 10;
ratio_array = NaN(length(cw_axis)-1,n_keep);
ir_array = NaN(length(cw_axis)-1,n_keep);
for c = 1:length(cw_axis)-1
    cw_filter = find(cw_vec_sharp>=cw_axis(c) & cw_vec_sharp<cw_axis(c+1) & sharpness_vec_sharp >=0);
    
    [~,si] = sort(ir_vec_sharp(cw_filter),'descend');
    
    rate_array = rate_array_sharp(cw_filter(si(1:n_keep)),:);
    
    ratio_array(c,:) = rate_array(:,9)./(rate_array(:,8)+rate_array(:,9)+rate_array(:,5));
    
    
    ir_array(c,:) = ir_vec_sharp(cw_filter(si(1:n_keep)));
end

x_array = repmat(1:length(cw_axis)-1,n_keep,1)';

%% identify maximally sharp and specific systems
%%%%%%%%%%
rate_bounds = [0.49 0.51];
prec_min = 0.95/sqrt(2);



metric_array_spec = sim_struct_spec.metric_array;   
rate_array_spec = sim_struct_spec.rate_array;   
sharpness_vec_spec = metric_array_spec(:,sharpness_index);
precision_vec_spec = sqrt(exp(metric_array_spec(:,precision_index)));
spec_vec_spec = metric_array_spec(:,spec_index);
pd_vec_spec = metric_array_spec(:,rate_index);


% identify candidate sharp system
sharp_options = find(pd_vec_sharp>=rate_bounds(1) & pd_vec_sharp<=rate_bounds(2) ...
                   & ir_vec_sharp >= prctile(ir_vec_sharp,99) & sharpness_vec_sharp>=0);
                 
[~, sharp_i] = sort(ir_vec_sharp(sharp_options),'descend');

sharp_rates = rate_array_sharp(sharp_options(sharp_i(1:end)),:);

% generate "swapped" version
swap_indices = [1 2 3 7 6 5 4 11 10 9 8];
sharp_rates_swap = sharp_rates(:,swap_indices);
metric_vec = calculateMetricsSym_v2(sharp_rates_swap, sim_info_sharp);                            












