% script to call core parameter sweep function to examine tradeoffs between
% different network behaviors

clear 
close all
addpath(genpath('../utilities/'))


[~,metric_names] = calculateMetricsNumeric([]);
nStates = 18;
% numerical_precision = 5;
functionPath = ['../utilities/metricFunctions/n' num2str(nStates) '_OR_NUM/'];

% make sure we're linked to the appropriate function subfolder
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));

% define save path
OutPath = ['../../out/bivariate_parameter_sweeps_n' num2str(nStates) '_numeric' filesep];
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
decision_rate_index = find(strcmp(metric_names,'DecisionRateNorm'));
decision_time_index = find(strcmp(metric_names,'DecisionTimeNorm'));
phi_index = find(strcmp(metric_names,'Phi'));
affinity_index = find(strcmp(metric_names,'AffinityVec'));
dev_index = find(strcmp(metric_names,'deviationFactor'));
cw_index = find(strcmp(metric_names,'CW'));


% set sim options
sweep_options = {'n_sim',1,'n_seeds',5,'n_iters_max',50,'nStates',nStates,'numCalcFlag',1,'numerical_precision',15};
%%
tic
[sim_info_eq, sim_struct_eq] = param_sweep_multi_v2([rate_index sharpness_index],functionPath, sweep_options{:},...
                                          'half_max_flag',false,'equilibrium_flag',false,'simType','all_specific');
toc;                                        
% [sim_info_neq, sim_struct_neq] = param_sweep_multi_v2([affinity_index spec_index],functionPath, sweep_options{:},...
%                                           'half_max_flag',false,'equilibrium_flag',false);                                        
%                                         
% toc                                        
%% %%%%%%%%%%%%%%%% info rate vs energy flux per cycle %%%%%%%%%%%%%%%%%%%%
tic
[sim_info_eq, sim_struct_eq] = param_sweep_multi_v2([cw_index decision_rate_index],functionPath, sweep_options{:},...
                                          'half_max_flag',false,'equilibrium_flag',true);
                                        
[sim_info_neq, sim_struct_neq] = param_sweep_multi_v2([cw_index decision_rate_index],functionPath, sweep_options{:},...
                                          'half_max_flag',false,'equilibrium_flag',false);                                        
                                                                                
toc     
