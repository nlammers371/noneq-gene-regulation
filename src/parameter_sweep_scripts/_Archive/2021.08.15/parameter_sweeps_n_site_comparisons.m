% script to call core parameter sweep function to examine tradeoffs between
% different network behaviors

clear 
close all
addpath(genpath('../utilities/'))

currentPath = pwd;
if strcmp(currentPath(1:7),'P:\Nick')
    DropboxPath = 'S:\Nick\Dropbox\Nonequilibrium\Nick\SweepOutput\';
else    
    DropboxPath = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\SweepOutput\';
end    
    
[~,metric_names] = calculateMetricsNumeric([]);
nStateVec = [6 18 54 162];
% numerical_precision = 5;
for n = 1:length(nStateVec)
    functionPathCell(n) = {['../utilities/metricFunctions/n' num2str(nStateVec(n)) '_OR_NUM/']};    
end

OutPath = [DropboxPath 'parameter_sweeps_multi_right_site\'];
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
ir_index = find(strcmp(metric_names,'DecisionRateNorm'));
decision_time_index = find(strcmp(metric_names,'DecisionTimeNorm'));
phi_index = find(strcmp(metric_names,'Phi'));
affinity_index = find(strcmp(metric_names,'AffinityVec'));
dev_index = find(strcmp(metric_names,'deviationFactor'));
cw_index = find(strcmp(metric_names,'CW'));


% set sim options
sweep_options = {'n_sim',50,'n_seeds',10,'n_iters_max',50,'numCalcFlag',1,'numerical_precision',10,'useParpool',1,'NumWorkers',8};

% IR vs Phi
metric_indices_ir = [phi_index ir_index];
results_struct_ir = struct;
for n = 1:length(nStateVec)-1
    tic
    [results_struct_ir(n).sim_info, results_struct_ir(n).sim_struct] = param_sweep_multi_v2(metric_indices_ir,functionPathCell{n}, sweep_options{:},...
                                              'half_max_flag', true,'equilibrium_flag', true,'simType','all_specific','nStates',nStateVec(n));
    results_struct_ir(n).functionPath =  functionPathCell{n};
    toc;    
end    

% save output
saveName = [metric_names{metric_indices_ir(1)} '_vs_' metric_names{metric_indices_ir(2)}];
save([OutPath saveName '_eq.mat'],'results_struct_ir', '-v7.3')

% Sharpness vs Precision
metric_indices_sh = [sharpness_index precision_index];
results_struct_sh = struct;
for n = 1:length(nStateVec)-1
    tic
    [results_struct_sh(n).sim_info, results_struct_sh(n).sim_struct] = param_sweep_multi_v2(metric_indices_sh,functionPathCell{n}, sweep_options{:},...
                                              'half_max_flag',true,'equilibrium_flag',true,'simType','all_specific','nStates',nStateVec(n));
    results_struct_sh(n).functionPath =  functionPathCell{n};
    toc;    
end    

% save output
saveName = [metric_names{metric_indices_sh(1)} '_vs_' metric_names{metric_indices_sh(2)}];
save([OutPath saveName '_eq.mat'],'results_struct_sh', '-v7.3')