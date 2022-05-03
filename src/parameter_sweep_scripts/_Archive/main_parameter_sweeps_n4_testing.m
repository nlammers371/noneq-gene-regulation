% script to call core parameter sweep function to examine tradeoffs between
% different network behaviors

clear 
close all
addpath(genpath('../utilities/'))

% set basic parameters
nStates = 4;
rate_bounds = repmat([-4 ; 4],1,3*nStates-4); % constrain transition rate magnitude
[~,metric_names] = calculateMetricsMultiState([]);

% make sure we're linked to the appropriate function subfolder
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(['../utilities/metricFunctions/n' num2str(nStates) '_OR_test/']));

% define save path
OutPath = ['../../out/bivariate_parameter_sweeps_n' num2str(nStates) '_OR_test' filesep];
mkdir(OutPath);
                         
% get index of useful metrics
flux_index = find(strcmp(metric_names,'Flux'));
phi_index = find(strcmp(metric_names,'Phi'));
rate_index = find(strcmp(metric_names,'Production Rate'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
precision_index = find(strcmp(metric_names,'Precision'));
decision_rate_index = find(strcmp(metric_names,'DecisionRateNorm'));
decision_time_index = find(strcmp(metric_names,'DecisionTimeNorm'));

% set sim options
sweep_options = {'n_seeds',5,'n_iters_max',50};
%% %%%%%%%%%%%%%%%%%%%% sharpness vs precision %%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
[sim_info_neq, sim_struct_neq] = param_sweep_multi([rate_index sharpness_index],sweep_options{:},...
                                                            'half_max_flag',false,'equilibrium_flag',false,'numTesting',1);
                                                          
% [sim_info_eq, sim_struct_eq] = param_sweep_multi([sharpness_index precision_index],sweep_options{:},...
%                                                             'half_max_flag',true,'equilibrium_flag',true);                                                          
                                                          
toc                                                          
%%
% set save name
save_name1 = ['param_sweep_results_' metric_names{sharpness_index} '_' ...
    metric_names{precision_index}];
save([OutPath save_name1 '_eq0.mat'],'sim_struct_neq','-v7.3')
save([OutPath save_name1 'info_eq0.mat'],'sim_info_neq')

save([OutPath save_name1 '_eq1.mat'],'sim_struct_eq','-v7.3')
save([OutPath save_name1 '_info_eq1.mat'],'sim_info_eq')

%% %%%%%%%%%%%%%%%% decision rate vs flux  %%%%%%%%%
tic
[sim_info_neq, sim_struct_neq]  = param_sweep_multi([flux_index decision_rate_index],sweep_options{:},...
              'half_max_flag',true,'equilibrium_flag',false);

toc

% save
save_name1 = ['param_sweep_results_' metric_names{phi_index} '_' ...
    metric_names{decision_rate_index}];
save([OutPath save_name1 '_eq0.mat'],'sim_struct_neq','-v7.3')
save([OutPath save_name1 'info_eq0.mat'],'sim_info_neq')

%% %%%%%%%%%%%%%%%%sharpness and precision vs flux %%%%%%%%%
tic
[sim_info_sharpness, sim_struct_sharpness]  = param_sweep_multi([flux_index sharpness_index],sweep_options{:},...
              'half_max_flag',true,'equilibrium_flag',false);

[sim_info_precision, sim_struct_precision]  = param_sweep_multi([flux_index precision_index],sweep_options{:},...
              'half_max_flag',true,'equilibrium_flag',false);            

toc

% save
save_name_sharp = ['param_sweep_results_' metric_names{flux_index} '_' ...
    metric_names{sharpness_index}];
save_name_precision = ['param_sweep_results_' metric_names{flux_index} '_' ...
    metric_names{precision_index}];
save([OutPath save_name_sharp '_eq0.mat'],'sim_struct_sharpness','-v7.3')
save([OutPath save_name_precision '_eq0.mat'],'sim_struct_precision')
