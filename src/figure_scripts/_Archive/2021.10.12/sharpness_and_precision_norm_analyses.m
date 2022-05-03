% script to examine validity of employing normalized definitions of
% sharpness and precision. Basically, I want to confirm that these
% proposed formulas do not depend on the transcription rate:
%   (a) in- and out-of equilibrium
%   (b) with and without Poisson noise from Pol II initiation

clear 
close all
addpath(genpath('../utilities/'))

% set basic parameters
nStates = 4;
rate_bounds = repmat([-8 ; 4],1,9); % constrain transition rate magnitude
[~,~,metric_names] = calculateMetricsSym_v2([]);

% make sure we're linked to the appropriate function subfolder
functionPath = ['../utilities/metricFunctions/n' num2str(nStates) '_OR/'];
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));

% define save path
OutPath = ['../../out/bivariate_parameter_sweeps_n' num2str(nStates) filesep];
mkdir(OutPath);
                         

% get index of useful metrics
rate_index = find(strcmp(metric_names,'Production Rate'));
sharpness_norm_index = find(strcmp(metric_names,'SharpnessNormed'));
precision_index = find(strcmp(metric_names,'Precision'));
precision_index_norm = find(strcmp(metric_names,'PrecisionNorm'));
precision_poisson_index = find(strcmp(metric_names,'PrecisionPoisson'));
pp_normed_index = find(strcmp(metric_names,'PrecPoissonNormed'));
tau_index = find(strcmp(metric_names,'CycleTime'));



% set sim options
sweep_options = {'n_seeds',10,'n_iters_max',50,'nStates',nStates,'numCalcFlag',0};

%% %%%%%%%%%%%%%%%% normalized sharpness vs precision %%%%%%%%%%%%%%%%%%%%
r_target = 1/3;
TauCycleLimit = 100;
tic
                                      
[simInfoNeq, sim_struct_neq] = param_sweep_multi_v3([sharpness_norm_index precision_index_norm],functionPath,sweep_options{:},...
                                          'half_max_flag',false,'equilibrium_flag',false,'r_target',r_target,'TauCycleLimit',TauCycleLimit);                                              
                                        
[simInfoEq, sim_struct_eq] = param_sweep_multi_v3([sharpness_norm_index precision_index_norm],functionPath,sweep_options{:},...
                                          'half_max_flag',false,'equilibrium_flag',true,'r_target',r_target,'TauCycleLimit',TauCycleLimit);                                              
toc     
%%
sharp_fig = figure;
hold on
scatter(sim_struct_neq(1).metric_array(:,sharpness_norm_index),exp(sim_struct_neq(1).metric_array(:,precision_index_norm)));
scatter(sim_struct_eq(1).metric_array(:,sharpness_norm_index),exp(sim_struct_eq(1).metric_array(:,precision_index_norm)));

%% %%%%%%%%%%%%%%%% normalized sharpness vs R %%%%%%%%%%%%%%%%%%%%
r_target = 1/3;
TauCycleLimit = 100;
tic
                                      
[simInfoNeq, sim_struct_neq] = param_sweep_multi_v3([rate_index sharpness_norm_index],functionPath,sweep_options{:},...
                                          'half_max_flag',false,'equilibrium_flag',false,'r_target',r_target,'TauCycleLimit',TauCycleLimit);                                              
                                        
[simInfoEq, sim_struct_eq] = param_sweep_multi_v3([rate_index sharpness_norm_index],functionPath,sweep_options{:},...
                                          'half_max_flag',false,'equilibrium_flag',true,'r_target',r_target,'TauCycleLimit',TauCycleLimit);                                              
toc     

sharp_fig = figure;
hold on
scatter(sim_struct_neq(1).metric_array(:,rate_index),sim_struct_neq(1).metric_array(:,sharp_norm_index));
scatter(sim_struct_eq(1).metric_array(:,rate_index),sim_struct_eq(1).metric_array(:,sharp_norm_index));


%% %%%%%%%%%%%%%%%% normalized precision vs R %%%%%%%%%%%%%%%%%%%%
close all
tic
                                      
[simInfoNeq, sim_struct_neq] = param_sweep_multi_v3([rate_index precision_index_norm],functionPath,sweep_options{:},...
                                          'half_max_flag',false,'equilibrium_flag',false,'r_target',r_target,'TauCycleLimit',TauCycleLimit);                                              
                                        
[simInfoEq, sim_struct_eq] = param_sweep_multi_v3([rate_index precision_index_norm],functionPath,sweep_options{:},...
                                          'half_max_flag',false,'equilibrium_flag',true,'r_target',r_target,'TauCycleLimit',TauCycleLimit);                                              
toc     

precision_fig = figure;
hold on
scatter(sim_struct_neq(1).metric_array(:,rate_index),exp(sim_struct_neq(1).metric_array(:,precision_index_norm)));
scatter(sim_struct_eq(1).metric_array(:,rate_index),exp(sim_struct_eq(1).metric_array(:,precision_index_norm)));


%% %%%%%%%%%%%%%%%% normalized precision (with Poisson noise) vs R %%%%%%%%%%%%%%%%%%%%
tic
                                      
[simInfoNeq, sim_struct_neq] = param_sweep_multi_v3([rate_index pp_norm_index],functionPath,sweep_options{:},...
                                          'half_max_flag',false,'equilibrium_flag',false,'r_target',r_target,'TauCycleLimit',TauCycleLimit);                                              
                                        
[simInfoEq, sim_struct_eq] = param_sweep_multi_v3([rate_index pp_norm_index],functionPath,sweep_options{:},...
                                          'half_max_flag',false,'equilibrium_flag',true,'r_target',r_target,'TauCycleLimit',TauCycleLimit);                                              
toc     

precision_fig = figure;
hold on
scatter(sim_struct_neq(1).metric_array(:,rate_index),exp(sim_struct_neq(1).metric_array(:,pp_norm_index)));
scatter(sim_struct_eq(1).metric_array(:,rate_index),exp(sim_struct_eq(1).metric_array(:,pp_norm_index)));
