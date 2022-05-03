% Plot results for IR vs energy for higher order models
clear 
close all
addpath(genpath('../utilities/'))

% %%%%%%%%%%%%%%%%  set relative read and write paths %%%%%%%%%%%%%%%%%%%%

DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\';
if ~exist(DropboxFolder)
    DropboxFolder = 'S:\Nick\Dropbox (Personal)\Nonequilibrium\Nick\';
end    

DataPath = [DropboxFolder  'SweepOutput\sweeps02_sharpness_vs_precision' filesep ];

% %%%%%%%%%%%%%%%%  Set plot parameters and constants %%%%%%%%%%%%%%%%%%%%
n_plot = 3e3; % number of points to plot
markerAlpha = 0.5; % marker transparency
markerSize = 75; % marker size

% get list of sweep results files with only 1 genera TF reaction
multi_bs_sweep_files = dir([DataPath 'sweep_results_s03_ns00_g01_cw0*']);
multi_bs_info_files = dir([DataPath 'sweep_info_s03_ns00_g01_cw0*']);

% load  
load([DataPath multi_bs_sweep_files(1).name])
load([DataPath multi_bs_info_files(1).name])

%%
[~,~,metric_names] = calculateMetricsNumeric_v3([]);
rate_index = find(strcmp(metric_names,'ProductionRate'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
precision_index = find(strcmp(metric_names,'Precision'));
phi_index = find(strcmp(metric_names,'Phi'));    
tau_index = find(strcmp(metric_names,'TauCycle'));    
ir_index = find(strcmp(metric_names,'IR'));


phi_filter = sim_results.metric_array(:,phi_index)>5;

% find sharp systems with r \approx 0.5
rate_filter = sim_results.metric_array(:,rate_index)>.4 & sim_results.metric_array(:,rate_index) <0.6;

options = find(rate_filter&phi_filter);

[sharp_max,sharp_i] = nanmax(sim_results.metric_array(options,precision_index))



%% calculate stuff
functionPath = sim_info.functionPath;
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));
numerical_precision = 10;

param_array = sim_results.rate_array(options(sharp_i),:);
n = 1;
valCellCS = mat2cell(param_array(n,:),size(param_array(n,:),1),ones(1,size(param_array(n,:),2)));    

Q_num = RSymFun(valCellCS{:});

% state probabilities
ss_vec = calculate_ss_num(Q_num,numerical_precision)
   