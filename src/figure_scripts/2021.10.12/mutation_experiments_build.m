% script to track information rate as a function of cw
clear 
close all
addpath(genpath('../utilities/'))

% set basic parameters
nStates = 6;
paramBounds = repmat([-10 ; 6],1,11); % constrain transition rate magnitude
[~,~,metric_names] = calculateMetricsSym_v2([]);

% specify function path
functionPath = ['../utilities/metricFunctions/n' num2str(nStates) '_OR/'];

% make sure we're linked to the appropriate function subfolder% make sure we're linked to the appropriate function subfolder
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));

% define save path

DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\manuscript\';
% DropboxFolder = 'S:\Nick\Dropbox\Nonequilibrium\Nick\manuscript\';

FigPath = [DropboxFolder 'observed_effects' filesep];
mkdir(FigPath);         

% get index of useful metrics
spec_index = find(strcmp(metric_names,'Specificity'));
sharp_right_norm_index = find(strcmp(metric_names,'SharpnessRightNorm'));
ir_index = find(strcmp(metric_names,'DecisionRateNorm'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
sharpness_norm_index = find(strcmp(metric_names,'SharpnessNormed'));
precision_index = find(strcmp(metric_names,'Precision'));
rate_index = find(strcmp(metric_names,'Production Rate'));
precision_right_index = find(strcmp(metric_names,'PrecisionRight'));
cw_index = find(strcmp(metric_names,'CW'));

% set sim options
sweep_options = {'n_seeds',5,'n_iters_max',50,'n_sim',10,'nStates',nStates};

% calculate sensitivity bound
alpha_factor = 100;

% specify plot options 
markerAlpha = 0.5; % marker transparency
markerSize = 75; % marker size
n_plot = 3e3;
bit_factor = log2(exp(1));

%% Run initial sweeps
cw_vec = 0:0.5:8;

% cw = 1e5; %
master_struct = struct;
for c = 1:length(cw_vec)
    tic
    [sim_info_neq, sim_struct_neq] = param_sweep_multi_v3([ir_index sharpness_index],functionPath,sweep_options{:},...
                                              'half_max_flag',true,'cw',10.^cw_vec(c),...
                                              'equilibrium_flag',false,'specFactor',alpha_factor,'paramBounds',paramBounds);
                                            
    [sim_info_eq, sim_struct_eq] = param_sweep_multi_v3([ir_index sharpness_index],functionPath,sweep_options{:},...
                                              'half_max_flag',true,'cw',10.^cw_vec(c),...
                                              'equilibrium_flag',true,'specFactor',alpha_factor,'paramBounds',paramBounds);  
                                            
    toc
                                            
    % extract key pieces of info
    
    % eq
    metric_array_eq = vertcat(sim_struct_eq.metric_array);
    rate_array_eq = vertcat(sim_struct_eq.rate_array);
    
    sharp_vec_eq = metric_array_eq(:,sharpness_index);
    ir_vec_eq = metric_array_eq(:,ir_index);
    ir_99_eq = prctile(ir_vec_eq,99.5);
    sharp_99_eq = prctile(sharp_vec_eq,99);
    optimal_ids_eq = find(ir_vec_eq>=ir_99_eq & sharp_vec_eq>=0);%>=sharp_99_eq);

%     master_struct(c).metric_array_eq = metric_array_eq;
%     master_struct(c).rate_array_eq = rate_array_eq;
    master_struct(c).sim_info_eq = sim_info_eq;
    master_struct(c).opt_rates_eq = rate_array_eq(optimal_ids_eq,:);
    master_struct(c).opt_ir_eq = ir_vec_eq(optimal_ids_eq);
    master_struct(c).opt_s_eq = sharp_vec_eq(optimal_ids_eq);
    master_struct(c).opt_ids_eq = optimal_ids_eq;
    
    master_struct(c).cw_val_eq = repmat(cw_vec(c),length(optimal_ids_eq),1);
    
    % neq (global)
    metric_array_neq = vertcat(sim_struct_neq.metric_array);
    rate_array_neq = vertcat(sim_struct_neq.rate_array);
    
    
    sharp_vec_neq = metric_array_neq(:,sharpness_index);
    ir_vec_neq = metric_array_neq(:,ir_index);
    ir_99_neq = prctile(ir_vec_neq,99.5);
    sharp_99_neq = prctile(sharp_vec_neq,99);
    optimal_ids_neq = find(ir_vec_neq>=ir_99_neq & sharp_vec_neq>=0);

%     master_struct(c).metric_array_neq = metric_array_neq;
    master_struct(c).sim_info_neq = sim_info_neq;
    master_struct(c).opt_rates_neq = rate_array_neq(optimal_ids_neq,:);
    master_struct(c).opt_ir_neq = ir_vec_neq(optimal_ids_neq);
    master_struct(c).opt_s_neq = sharp_vec_neq(optimal_ids_neq);
    master_struct(c).opt_ids_neq = optimal_ids_neq;
                    
    % store cw values
    master_struct(c).cw_val_neq = repmat(cw_vec(c),length(optimal_ids_neq),1);
end

save([FigPath 'mutation_experiment_data'],'master_struct')