% script to generate data for system sharpness behaviors
clear 
close all
addpath(genpath('../utilities/'))
% define save path
OutPath = ['../../out/bivariate_parameter_sweeps/'];
mkdir(OutPath);

rate_bounds = [-4 -4 -4 -4 ;...
                4  4  4  4 ];%repmat([-4 ; 1],1,8); % constrain transition rate magnitude
              
rate_bounds = repmat(rate_bounds,1,2);              
[~,metric_names] = calculateMetrics_v4([]);

% get index of useful metrics
flux_index = find(strcmp(metric_names,'Flux'));
rate_index = find(strcmp(metric_names,'Production Rate'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
sharpness_norm_index = find(strcmp(metric_names,'SharpnessNormed'));
precision_index = find(strcmp(metric_names,'Precision'));
decision_rate_index = find(strcmp(metric_names,'DecisionRateNorm'));
decision_time_index = find(strcmp(metric_names,'DecisionTimeNorm'));
affinity_index = find(strcmp(metric_names,'Affinity'));
switch_noise_index = find(strcmp(metric_names,'SwitchingNoise'));
bi_noise_index = find(strcmp(metric_names,'BinomialNoise'));
neq_noise_factor_index = find(strcmp(metric_names,'NeqNoiseFactor'));

% set sim options
sweep_options = {'n_seeds',5,'n_iters_max',50,'rate_bounds',rate_bounds};

%% %%%%%%%%%%%%%%%% info rate vs energy flux per cycle %%%%%%%%%%%%%%%%%%%%
tic
[sim_info_neq, sim_struct_neq] = param_sweep_v5([flux_index decision_rate_index],sweep_options{:},'half_max_flag',false);

[sim_info_neq_half, sim_struct_neq_half] = param_sweep_v5([flux_index decision_rate_index],sweep_options{:},'half_max_flag',true);
toc

% set save name
save_name_flux = ['param_sweep_results_' metric_names{flux_index} '_' ...
    metric_names{decision_rate_index}];
save([OutPath save_name_flux '_eq0.mat'],'sim_struct_neq','-v7.3')
save([OutPath save_name_flux 'info_eq0.mat'],'sim_info_neq')

save([OutPath save_name_flux '_half_eq0.mat'],'sim_struct_neq','-v7.3')
save([OutPath save_name_flux '_half_info_eq0.mat'],'sim_info_neq')

%% %%%%%%%%%%%%%%%% info rate vs "affinity" (i.e. log(koff/kon) ) %%%%%%%%%
tic
[sim_info_neq, sim_struct_neq]  = param_sweep_v5([affinity_index decision_rate_index],sweep_options{:},'equilibrium_flag',false);

[sim_info_eq, sim_struct_eq]  = param_sweep_v5([affinity_index decision_rate_index],sweep_options{:},'equilibrium_flag',true);
% sim_struct_neq = param_sweep_v5([affinity_index decision_rate_index],sweep_options{:},'equilibrium_flag',false,'n_seeds',50);
toc

% set save name
save_name_affinity = ['param_sweep_results_' metric_names{affinity_index} '_' ...
    metric_names{decision_rate_index}];
save([OutPath save_name_affinity '_eq0.mat'],'sim_struct_neq','-v7.3')
save([OutPath save_name_affinity 'info_eq0.mat'],'sim_info_neq')
save([OutPath save_name_affinity '_eq1.mat'],'sim_struct_eq','-v7.3')
save([OutPath save_name_affinity 'info_eq1.mat'],'sim_info_eq')

%% %%%%%%%%%%%%%%%% decision time vs "affinity" (i.e. log(koff/kon) ) %%%%%%%%%
tic
[sim_info_neq, sim_struct_neq] = param_sweep_v5([affinity_index decision_time_index],sweep_options{:},'equilibrium_flag',false);

[sim_info_eq, sim_struct_eq] = param_sweep_v5([affinity_index decision_time_index],sweep_options{:},'equilibrium_flag',true);
% sim_struct_neq = param_sweep_v5([affinity_index decision_rate_index],sweep_options{:},'equilibrium_flag',false,'n_seeds',50);
toc

% set save name
save_name_affinity = ['param_sweep_results_' metric_names{affinity_index} '_' ...
    metric_names{decision_time_index}];
save([OutPath save_name_affinity '_eq0.mat'],'sim_struct_neq','-v7.3')
save([OutPath save_name_affinity 'info_eq0.mat'],'sim_info_neq')
save([OutPath save_name_affinity '_eq1.mat'],'sim_struct_eq','-v7.3')
save([OutPath save_name_affinity 'info_eq1.mat'],'sim_info_eq')


%% %%%%%%%%%%%%%%%% sharpness vs precision %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
[sim_info_neq, sim_struct_neq] = param_sweep_v5([sharpness_index precision_index],sweep_options{:},...
                        'equilibrium_flag',false,'half_max_flag',false);
toc

tic
[sim_info_eq, sim_struct_eq] = param_sweep_v5([sharpness_index precision_index],sweep_options{:},...
                        'equilibrium_flag',true,'half_max_flag',false);
toc

% set save name
save_name_tradeoff = ['param_sweep_results_' metric_names{sharpness_index} '_' ...
                      metric_names{precision_index}];
save([OutPath save_name_tradeoff '_eq0.mat'],'sim_struct_neq','-v7.3')
save([OutPath save_name_tradeoff 'info_eq0.mat'],'sim_info_neq')
save([OutPath save_name_tradeoff '_eq1.mat'],'sim_struct_eq','-v7.3')
save([OutPath save_name_tradeoff 'info_eq1.mat'],'sim_info_eq')

%% %%%%%%%%%%%%%%%% sharpness vs switching noise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_indices = [sharpness_index switch_noise_index];
tic
[sim_info_neq, sim_struct_neq] = param_sweep_v5(plot_indices, sweep_options{:},...
                        'equilibrium_flag',false,'half_max_flag',true);
toc

tic
[sim_info_eq, sim_struct_eq] = param_sweep_v5(plot_indices,sweep_options{:},...
                        'equilibrium_flag',true,'half_max_flag',true);
toc

% set save name
save_name_tradeoff = ['param_sweep_results_' metric_names{plot_indices(1)} '_' ...
                      metric_names{plot_indices(2)}];
save([OutPath save_name_tradeoff '_eq0.mat'],'sim_struct_neq','-v7.3')
save([OutPath save_name_tradeoff 'info_eq0.mat'],'sim_info_neq')
save([OutPath save_name_tradeoff '_eq1.mat'],'sim_struct_eq','-v7.3')
save([OutPath save_name_tradeoff 'info_eq1.mat'],'sim_info_eq')

%% %%%%%%%%%%% sharpness vs precision (with half-max constraint) %%%%%%%%%%
plot_indices = [sharpness_index switch_noise_index];

tic
[sim_info_neq, sim_struct_neq] = param_sweep_v5(plot_indices,sweep_options{:},...
                        'equilibrium_flag',false,'half_max_flag',true);
toc

tic
[sim_info_eq, sim_struct_eq] = param_sweep_v5(plot_indices,sweep_options{:},...
                        'equilibrium_flag',true,'half_max_flag',true);
toc

% set save name
save_name_tradeoff = ['param_sweep_results_' metric_names{plot_indices(1)} '_' ...
                      metric_names{plot_indices(2)}];
save([OutPath save_name_tradeoff '_half_eq0.mat'],'sim_struct_neq','-v7.3')
save([OutPath save_name_tradeoff '_halfinfo_eq0.mat'],'sim_info_neq')
save([OutPath save_name_tradeoff '_half_eq1.mat'],'sim_struct_eq','-v7.3')
save([OutPath save_name_tradeoff '_halfinfo_eq1.mat'],'sim_info_eq')

%% %%%%%%%%%%%%%%%% pon vs info rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_indices = [sharpness_norm_index neq_noise_factor_index];
tic
[sim_info_neq, sim_struct_neq] = param_sweep_v5(plot_indices, sweep_options{:},...
                        'equilibrium_flag',false,'half_max_flag',false);
toc

tic
[sim_info_eq, sim_struct_eq] = param_sweep_v5(plot_indices,sweep_options{:},...
                        'equilibrium_flag',true,'half_max_flag',false);
toc

% set save name
save_name_tradeoff = ['param_sweep_results_' metric_names{plot_indices(1)} '_' ...
                      metric_names{plot_indices(2)}];
save([OutPath save_name_tradeoff '_eq0.mat'],'sim_struct_neq','-v7.3')
save([OutPath save_name_tradeoff 'info_eq0.mat'],'sim_info_neq')
save([OutPath save_name_tradeoff '_eq1.mat'],'sim_struct_eq','-v7.3')
save([OutPath save_name_tradeoff 'info_eq1.mat'],'sim_info_eq')