% script to track observed sharpness as a function of CW
clear 
close all
addpath(genpath('../utilities/'))

% set basic parameters
nStates = 4;
paramBounds = repmat([-10 ; 6],1,9); % constrain transition rate magnitude
[~,~,metric_names] = calculateMetricsSym_v2([]);

% specify function path
functionPath = ['../utilities/metricFunctions/n' num2str(nStates) '_OR/'];

% make sure we're linked to the appropriate function subfolder% make sure we're linked to the appropriate function subfolder
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));

% define save path

DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\manuscript\';
% DropboxFolder = 'S:\Nick\Dropbox\Nonequilibrium\Nick\manuscript\';

FigPath = [DropboxFolder 'exploratory_analyses' filesep];
mkdir(FigPath);         

% get index of useful metrics
spec_index = find(strcmp(metric_names,'Specificity'));
sharp_right_norm_index = find(strcmp(metric_names,'SharpnessRightNorm'));
ir_index = find(strcmp(metric_names,'DecisionRateNorm'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
sharpness_norm_index = find(strcmp(metric_names,'SharpnessNormed'));
precision_index = find(strcmp(metric_names,'Precision'));
precision_norm_index = find(strcmp(metric_names,'PrecisionNorm'));
rate_index = find(strcmp(metric_names,'Production Rate'));
precision_right_index = find(strcmp(metric_names,'PrecisionRight'));
cw_index = find(strcmp(metric_names,'CW'));
precision_poisson_index = find(strcmp(metric_names,'PrecisionPoisson'));
ir_poisson_index = find(strcmp(metric_names,'DecisionRatePoisson'));
phi_index = find(strcmp(metric_names,'Phi'));

% set sim options
sweep_options = {'n_seeds',5,'n_iters_max',50,'n_sim',5,'nStates',nStates};

% calculate sensitivity bound
alpha_factor = 100;
% seq = 1/4;

% specify plot options 
markerAlpha = 0.5; % marker transparency
markerSize = 75; % marker size
% n_plot = 3e3;
% bit_factor = log2(exp(1));
%% 
% observed sharpness vs cw
tic
[sim_info_neq, sim_struct_neq] = param_sweep_multi_v3([sharpness_index precision_poisson_index],functionPath,sweep_options{:},...
                                          'half_max_flag',false,...
                                          'equilibrium_flag',false,'paramBounds',paramBounds,'TauCycleLimit',5*60);
                                        
[sim_info_eq, sim_struct_eq] = param_sweep_multi_v3([sharpness_index precision_poisson_index],functionPath,sweep_options{:},...
                                          'half_max_flag',false,...
                                          'equilibrium_flag',true,'paramBounds',paramBounds,'TauCycleLimit',5*60);     
                                        
[sim_info_ir, sim_struct_ir] = param_sweep_multi_v3([phi_index ir_poisson_index],functionPath,sweep_options{:},...
                                          'half_max_flag',false,...
                                          'equilibrium_flag',false,'paramBounds',paramBounds,'TauCycleLimit',5*60);                                        
                                       
toc     

%% Let's pull out a random set of rates and see what happens when we fiddle with stuff
% close all
metric_array_neq = vertcat(sim_struct_neq.metric_array);
pp_vec = exp(metric_array_neq(:,precision_poisson_index));
s_vec = metric_array_neq(:,sharpness_index);
ir_vec = metric_array_neq(:,ir_poisson_index);
r_vec = metric_array_neq(:,rate_index);

metric_array_eq = vertcat(sim_struct_eq.metric_array);
pp_vec_eq = exp(metric_array_eq(:,precision_poisson_index));
s_vec_eq = metric_array_eq(:,sharpness_index);
ir_vec_eq = metric_array_eq(:,ir_poisson_index);
r_vec_eq = metric_array_eq(:,rate_index);

figure;
hold on
scatter(s_vec./(r_vec.*(1-r_vec)),pp_vec.*(r_vec.*(1-r_vec)).^2,[],ir_vec)
scatter(s_vec_eq./(r_vec_eq.*(1-r_vec_eq)),pp_vec_eq.*(r_vec_eq.*(1-r_vec_eq)).^2)
% set(gca,'yscale','log')
% ylim([0 20])
colorbar

%%
metric_array_ir = vertcat(sim_struct_ir.metric_array);
phi_vec = metric_array_ir(:,phi_index);
ir_vec_ir = metric_array_ir(:,ir_poisson_index);

figure;
hold on
scatter(phi_vec,ir_vec_ir)
set(gca,'xscale','log')
xlim([1e-5 5e1])

%% Compare bounds implied by normalized sweep to true HM results

tic
[sim_info_hm, sim_struct_hm] = param_sweep_multi_v3([sharpness_index precision_index],functionPath,sweep_options{:},...
                                          'half_max_flag',true,...
                                          'equilibrium_flag',false,'paramBounds',paramBounds,'TauCycleLimit',5*60);
                                        
[sim_info_norm, sim_struct_norm] = param_sweep_multi_v3([sharpness_norm_index precision_norm_index],functionPath,sweep_options{:},...
                                          'half_max_flag',false,...
                                          'equilibrium_flag',false,'paramBounds',paramBounds,'TauCycleLimit',5*60);
                                        
                                        
%% It looks like thelimits are accurate, even if the implied metric values for sub-optimal networks are not
close all
figure;
hold on
scatter(sim_struct_hm(1).metric_array(:,sharpness_index),exp(sim_struct_hm(1).metric_array(:,precision_index)))
scatter(sim_struct_norm(1).metric_array(:,sharpness_norm_index)/4,exp(sim_struct_norm(1).metric_array(:,precision_norm_index))*16)

xlabel('sharpness')
ylabel('precision')

%% Confirm that IR vs. Phi curve is HM-independent

[sim_info_hm_ir, sim_struct_hm_ir] = param_sweep_multi_v3([phi_index ir_index],functionPath,sweep_options{:},...
                                          'half_max_flag',true,...
                                          'equilibrium_flag',false,'paramBounds',paramBounds,'TauCycleLimit',5*60);
                                        
[sim_info_ir, sim_struct_ir] = param_sweep_multi_v3([phi_index ir_index],functionPath,sweep_options{:},...
                                          'half_max_flag',false,...
                                          'equilibrium_flag',false,'paramBounds',paramBounds,'TauCycleLimit',5*60);      
%%
close all
figure;
hold on
scatter(sim_struct_hm_ir(1).metric_array(:,phi_index),sim_struct_hm_ir(1).metric_array(:,ir_index))
scatter(sim_struct_ir(1).metric_array(:,phi_index),sim_struct_ir(1).metric_array(:,ir_index))

set(gca,'xscale','log')
xlim([1e-5 5e1])
xlabel('energy dissipation rate (\Phi)')
ylabel('information rate (IR)')