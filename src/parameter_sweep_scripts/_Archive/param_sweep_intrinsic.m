% script to generate data and figures for fidelity v sharpness parameter
% space exploration
clear 
close all
addpath(genpath('../utilities/'));
OutPath = ['../../out/bivariate_parameter_sweeps/'];

[~,metric_names] = calculateMetrics_v4([]);
flux_index = find(strcmp(metric_names,'Flux'));

% call edge sampler for production rate vs. sharpness
metric_indices_int = [1,4];
metric_indices_flux = [flux_index,5];

% constrain transition rate magnitud
rate_bounds = [-4 4];%log10([1e-4, 1e4]); 

options = {'rate_bounds',rate_bounds,'constrain_off_rates',false,'half_max_flag',true};
% run simulation for sharpness vs intrinsic noise (out of equiilibrium)
tic
sim_struct_int = param_sweep_v5(metric_indices_int,'n_seeds',10,'n_iters',50,...
  'rnd_seed',435,options{:});
toc

% run simulation for sharpness vs intrinsic noise (at equiilibrium)
tic
sim_struct_int_eq = param_sweep_v5(metric_indices_int,'n_seeds',5,'equilibrium_flag',1,...
  'rnd_seed',122,options{:});
toc

% save
save_name_int = ['param_sweep_results_' metric_names{metric_indices_int(1)} '_' ...
    metric_names{metric_indices_int(2)}];
save([OutPath save_name_int '_eq0.mat'],'sim_struct_int','-v7.3')
save([OutPath save_name_int '_eq1.mat'],'sim_struct_int_eq','-v7.3')

% "idealized" information vs flux
%%
tic
sim_struct_flux = param_sweep_v5(metric_indices_flux,'n_seeds',10,'n_iters',50,...
  'rnd_seed',444,options{:});
toc

% set save name
%%
save_name_flux = ['param_sweep_results_' metric_names{metric_indices_flux(1)} '_' ...
    metric_names{metric_indices_flux(2)} '_eq0.mat'];
save([OutPath save_name_flux],'sim_struct_flux','-v7.3')
