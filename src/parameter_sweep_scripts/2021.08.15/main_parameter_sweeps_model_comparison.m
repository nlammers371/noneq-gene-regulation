% script to call core parameter sweep function to examine tradeoffs between
% different network behaviors

clear 
close all
addpath(genpath('../utilities/'))

[~,metric_names_num] = calculateMetricsNumeric([]);
[~,metric_names_sym] = calculateMetricsSym([]);
nStateVec = [6];% 6 12];
folderRoot = '../utilities/metricFunctions/';
folderNameCell = {'n6_OR_Nucleosome_NUM'};
numCalcFlags = [1];% 0 1];
simTypeCell = {'Nucleosome', 'PIC', 'General'};
sweepVarIndexCell = cell(1,length(nStateVec));
% numerical_precision = 5;

% get index of useful metrics
flux_index = find(strcmp(metric_names_num,'Flux'));
rate_index = find(strcmp(metric_names_num,'Production Rate'));
spec_index = find(strcmp(metric_names_num,'Specificity'));
spec_alt_index = find(strcmp(metric_names_num,'specFactorAlt'));
sharp_right_index = find(strcmp(metric_names_num,'SharpnessRight'));
sharpness_index = find(strcmp(metric_names_num,'Sharpness'));
sharpness_norm_index = find(strcmp(metric_names_num,'SharpnessNormed'));
precision_index = find(strcmp(metric_names_num,'Precision'));
ir_index = find(strcmp(metric_names_num,'DecisionRateNorm'));
decision_time_index = find(strcmp(metric_names_num,'DecisionTimeNorm'));
phi_index = find(strcmp(metric_names_num,'Phi'));
affinity_index = find(strcmp(metric_names_num,'AffinityVec'));
dev_index = find(strcmp(metric_names_num,'deviationFactor'));
cw_index = find(strcmp(metric_names_num,'CW'));

% set sim options
sweep_options = {'n_sim',1,'n_seeds',5,'n_iters_max',50,'numerical_precision',15,'numCalcFlag',1};
%%
results_struct = struct;
for n = 1:length(folderNameCell)
    functionPath  = [folderRoot folderNameCell{n}];
    tic
    [results_struct(n).sim_info, results_struct(n).sim_struct] = param_sweep_multi_v2([sharpness_index, precision_index],functionPath, ...
                                              sweep_options{:},'half_max_flag',false,'equilibrium_flag',false,...
                                              'nStates',nStateVec(n));
                                            
    
%     [results_struct(n).sim_info_hm, results_struct(n).sim_struct_hm] = param_sweep_multi_v2([phi_index,ir_index],functionPath, ...
%                                               sweep_options{:},'half_max_flag',true,'equilibrium_flag',false,...
%                                               'nStates',nStateVec(n));
    toc;    
end    
 
%%
% close all
figure;
hold on
r_vec = results_struct(n).sim_struct.metric_array(:,rate_index);
scatter(results_struct(n).sim_struct.metric_array(:,sharpness_index) ./ (r_vec.*(1-r_vec)),...
        exp(results_struct(n).sim_struct.metric_array(:,precision_index)).^2,[],results_struct(n).sim_struct.metric_array(:,ir_index))
% scatter(results_struct(n).sim_struct.metric_array(:,phi_index),...
%         results_struct(n).sim_struct.metric_array(:,ir_index))
      
% scatter(results_struct(n).sim_struct_hm.metric_array(:,phi_index),...
%         results_struct(n).sim_struct_hm.metric_array(:,ir_index))      