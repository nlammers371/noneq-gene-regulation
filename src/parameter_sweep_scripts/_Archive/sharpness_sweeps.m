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

flux_index = find(strcmp(metric_names,'Flux'));
rate_index = find(strcmp(metric_names,'Production Rate'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
sharpness_norm_index = find(strcmp(metric_names,'SharpnessNormed'));
% call edge sampler for production rate vs. sharpness

metric_indices_flux = [flux_index,1];


% set sim options
sweep_options = {'n_seeds',5,'n_iters_max',50, 'rate_bounds',rate_bounds};

%% %%%%%%%%%%%%%%%% run simulation for pd rate vs sharpness %%%%%%%%%%%%%%%
metric_indices_1 = [rate_index, sharpness_index];

tic
[sim_info_neq, sim_results_neq] = param_sweep_v5(metric_indices_1,sweep_options{:});
toc

% now at equilibrium
tic
[sim_info_eq, sim_results_eq] = param_sweep_v5(metric_indices_1,sweep_options{:},'equilibrium_flag',true);
toc

% set save name
save_name_pd = ['param_sweep_results_' metric_names{metric_indices_1(1)} '_' ...
    metric_names{metric_indices_1(2)}];
save([OutPath save_name_pd '_eq0.mat'],'sim_results_neq','-v7.3')
save([OutPath save_name_pd '_info_eq0.mat'],'sim_info_neq','-v7.3')
save([OutPath save_name_pd '_eq1.mat'],'sim_results_eq','-v7.3')
save([OutPath save_name_pd '_info_eq1.mat'],'sim_info_eq','-v7.3')

%% %%%%%%%%%%%%%%%%%%%%%%%%%% flux vs sharpness %%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
[sim_info, sim_results] = param_sweep_v5([flux_index sharpness_index],sweep_options{:});
toc

% now with half-max constraint
tic
[sim_info_half, sim_results_half] = param_sweep_v5([flux_index sharpness_index],sweep_options{:},'half_max_flag',true);
toc

% set save name
save_name_flux = ['param_sweep_results_' metric_names{flux_index} '_' ...
    metric_names{sharpness_index}];
save([OutPath save_name_flux '_eq0.mat'],'sim_results','-v7.3')
save([OutPath save_name_flux '_info_eq0.mat'],'sim_info','-v7.3')
save([OutPath save_name_flux '_half_eq0.mat'],'sim_results_half','-v7.3')
save([OutPath save_name_flux '_half_info_eq0.mat'],'sim_info_half','-v7.3')


%% %%%%%%%%%%%%%%%% sharpness vs fraction of time in 1 or 3 %%%%%%%%%%%%%%%

tic
sim_struct_frac_neq = param_sweep_v5([sharpness_norm_index frac_index],sweep_options{:},'equilibrium_flag',false);

sim_struct_frac_eq = param_sweep_v5([sharpness_norm_index frac_index],sweep_options{:},'equilibrium_flag',true);
toc

% set save name
save_name_frac = ['param_sweep_results_' metric_names{sharpness_norm_index} '_' ...
    metric_names{frac_index}];
save([OutPath save_name_frac '_eq0.mat'],'sim_struct_frac_neq','-v7.3')
save([OutPath save_name_frac '_eq1.mat'],'sim_struct_frac_eq','-v7.3')



%% %%%%%%%%%%%%%%%% dkon/dc vs dkoff/dc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dkon_index = find(strcmp(metric_names,'dKondC'));
dkoff_index = find(strcmp(metric_names,'dKoffdC'));
dkdc_indices = [dkon_index dkoff_index];
% set save name
save_name_dkdc = ['param_sweep_results_' metric_names{dkdc_indices(1)} '_' ...
    metric_names{dkdc_indices(2)}];

tic
sim_struct = param_sweep_v5(dkdc_indices,sweep_options{:});
save([OutPath save_name_dkdc '_eq0.mat'],'sim_struct','-v7.3')
toc

% now at equilibrium
tic
sim_struct = param_sweep_v5([dkon_index dkoff_index],sweep_options{:},'equilibrium_flag',true);
toc
save([OutPath save_name_dkdc '_eq1.mat'],'sim_struct','-v7.3')


%% sharpness vs flux
tic
sim_struct_flux = param_sweep_v5(metric_indices_flux,'n_seeds',5,'n_iters',50,...
  'half_max_flag',true,sweep_options{:});
toc
% set save name
save_name_flux = ['param_sweep_results_' metric_names{metric_indices_flux(1)} '_' ...
    metric_names{metric_indices_flux(2)} '_eq0.mat'];
save([OutPath save_name_flux],'sim_struct_flux','-v7.3')

% %%
% pdNeq = reshape(sim_struct_pd.metric_array(:,2,:),1,[]);
% shNeq = reshape(sim_struct_pd.metric_array(:,1,:),1,[]);
% koffNeq = reshape(sim_struct_pd.metric_array(:,14,:),1,[]);
% konNeq = reshape(sim_struct_pd.metric_array(:,13,:),1,[]);
% dKoffdCNeq = reshape(sim_struct_pd.metric_array(:,16,:),1,[]);
% dKondCNeq = reshape(sim_struct_pd.metric_array(:,15,:),1,[]);
% opt_indices_neq = find(shNeq >= 1.9*pdNeq.*(1-pdNeq));
% % opt_indices_neq = find(shNeq >= 0.47);
% 
% %%
% pdEq = reshape(sim_struct_pd_eq.metric_array(:,2,:),1,[]);
% shEq = reshape(sim_struct_pd_eq.metric_array(:,1,:),1,[]);
% koffEq = reshape(sim_struct_pd_eq.metric_array(:,14,:),1,[]);
% konEq = reshape(sim_struct_pd_eq.metric_array(:,13,:),1,[]);
% dKoffdCEq = reshape(sim_struct_pd_eq.metric_array(:,16,:),1,[]);
% dKondCEq = reshape(sim_struct_pd_eq.metric_array(:,15,:),1,[]);
% opt_indices_eq = find(shEq >= 0.95*pdEq.*(1-pdEq));
% % opt_indices_eq = find(shEq >= 0.23);