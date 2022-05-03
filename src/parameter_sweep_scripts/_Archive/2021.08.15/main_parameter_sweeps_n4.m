% script to call core parameter sweep function to examine tradeoffs between
% different network behaviors

clear 
close all
addpath(genpath('../utilities/'))


[~,metric_names] = calculateMetricsSym([]);
nStates = 4;
functionPath = ['../utilities/metricFunctions/n' num2str(nStates) '_OR/'];

% make sure we're linked to the appropriate function subfolder
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));

% define save path
OutPath = ['../../out/bivariate_parameter_sweeps_n' num2str(nStates) '_OR' filesep];
mkdir(OutPath);
                         

% get index of useful metrics
flux_index = find(strcmp(metric_names,'Flux'));
rate_index = find(strcmp(metric_names,'Production Rate'));
spec_index = find(strcmp(metric_names,'Specificity'));
spec_alt_index = find(strcmp(metric_names,'specFactorAlt'));
sharp_right_index = find(strcmp(metric_names,'SharpnessRight'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
sharpness_norm_index = find(strcmp(metric_names,'SharpnessNormed'));
precision_index = find(strcmp(metric_names,'Precision'));
information_rate_index = find(strcmp(metric_names,'DecisionRateNorm'));
decision_time_index = find(strcmp(metric_names,'DecisionTimeNorm'));
phi_index = find(strcmp(metric_names,'Phi'));
affinity_index = find(strcmp(metric_names,'AffinityVec'));
dev_index = find(strcmp(metric_names,'deviationFactor'));

% set sim options
sweep_options = {'n_seeds',5,'n_iters_max',50,'nStates',nStates,'numCalcFlag',0};

%% %%%%%%%%%%%%%%%%%%%% IR vs Flux %%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
[simInfoIR, simStructIR] = param_sweep_multi_v2([phi_index information_rate_index],functionPath,sweep_options{:},...
                                                            'half_max_flag',true,'equilibrium_flag',false);                                                                                                               
                                                          
toc                                                          
% set save name
save_name_ir = ['param_sweep_results_' metric_names{phi_index} '_' metric_names{information_rate_index}];

save([OutPath save_name_ir '_eq0.mat'],'simStructIR','-v7.3')
save([OutPath save_name_ir '_info_eq0.mat'],'simInfoIR')

%% %%%%%%%%%%%%%%%% decision rate vs flux  %%%%%%%%%
tic
[sim_info_neq, sim_struct_neq]  = param_sweep_multi([flux_index information_rate_index],sweep_options{:},...
              'half_max_flag',true,'equilibrium_flag',false);

toc

% save
save_name_ir = ['param_sweep_results_' metric_names{phi_index} '_' ...
    metric_names{information_rate_index}];
save([OutPath save_name_ir '_eq0.mat'],'sim_struct_neq','-v7.3')
save([OutPath save_name_ir 'info_eq0.mat'],'sim_info_neq')

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

%%

tic
[sim_info_sharpness, sim_struct_sharpness]  = param_sweep_multi([flux_index sharpness_index],sweep_options{:},...
              'half_max_flag',true,'equilibrium_flag',false);

[sim_info_precision, sim_struct_precision]  = param_sweep_multi([flux_index precision_index],sweep_options{:},...
              'half_max_flag',true,'equilibrium_flag',false);            

toc

