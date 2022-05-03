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
    
[~,~,metric_names] = calculateMetricsNumeric_v3([]);
nStateVec = [6 18 54 162];
% numerical_precision = 5;
for n = 1:length(nStateVec)
    functionPathCell(n) = {['../utilities/metricFunctions/n' num2str(nStateVec(n)) '_OR_NUM/']};    
end

OutPath = [DropboxPath 'parameter_sweeps_multi_wrong_site_test/'];
mkdir(OutPath);
                         

% get index of useful metrics
flux_index = find(strcmp(metric_names,'Flux'));
rate_index = find(strcmp(metric_names,'Production Rate'));
spec_index = find(strcmp(metric_names,'Specificity'));
spec_alt_index = find(strcmp(metric_names,'specFactorAlt'));
spec_full_index = find(strcmp(metric_names,'specificityFactorFull'));
sharp_right_index = find(strcmp(metric_names,'SharpnessRight'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
sharpness_norm_index = find(strcmp(metric_names,'SharpnessRightNorm'));
precision_index = find(strcmp(metric_names,'Precision'));
ir_index = find(strcmp(metric_names,'DecisionRateNorm'));
decision_time_index = find(strcmp(metric_names,'DecisionTimeNorm'));
phi_index = find(strcmp(metric_names,'Phi'));
affinity_index = find(strcmp(metric_names,'AffinityVec'));
dev_index = find(strcmp(metric_names,'deviationFactor'));
cw_index = find(strcmp(metric_names,'CW'));


% set sim options
sweep_options = {'n_sim',1,'n_seeds',10,'n_iters_max',50,'numCalcFlag',1,'numerical_precision',12,'useParpool',0,'NumWorkers',4};

%% IR vs CW for each motif
% metric_indices_ir = [cw_index ir_index];

wb = waitbar(0,'Running IR parameter sweeps for different motifs...');
totalRuns = length(nStateVec)*3;
for n = 1%length(nStateVec)
    results_struct_ir = struct;
%     tic
%     % sharpness-optimized
%     [results_struct_ir.sim_info_sh, results_struct_ir.sim_struct_sh] = param_sweep_multi_v2([cw_index sharpness_norm_index],...
%                                               functionPathCell{n}, sweep_options{:},...
%                                               'half_max_flag', false,'equilibrium_flag',false,'nStates',nStateVec(n));
%     results_struct_ir.sharp_opt_indices = [cw_index sharpness_norm_index];
%     toc
%     
%     waitbar(((n-1)*length(nStateVec)+1)/totalRuns,wb)
%   
%     tic
%     % specificity-optimized
%     [results_struct_ir.sim_info_sp, results_struct_ir.sim_struct_sp] = param_sweep_multi_v2([cw_index spec_index],...
%                                               functionPathCell{n}, sweep_options{:},...
%                                               'half_max_flag', false,'equilibrium_flag',false,'nStates',nStateVec(n));
%     results_struct_ir.spec_opt_indices = [cw_index spec_index];
%     toc
%     waitbar(((n-1)*length(nStateVec)+2)/totalRuns,wb)
%     tic
%     
    % IR-optimized (global optimum)
    [results_struct_ir.sim_info_ir, results_struct_ir.sim_struct_ir] = param_sweep_multi_v3([cw_index ir_index],...
                                              functionPathCell{n}, sweep_options{:},...
                                              'half_max_flag', false,'equilibrium_flag',false,'nStates',nStateVec(n));
    results_struct_ir.ir_opt_indices = [cw_index ir_index];
    
    results_struct_ir.functionPath =  functionPathCell{n};
    waitbar(((n-1)*length(nStateVec)+3)/totalRuns,wb)
%     
%     toc;    
    % save output
    disp(['Saving results for ' num2str(nStateVec(n)) ' states...'])
    save([OutPath 'ir_sharp_spec_n'  num2str(nStateVec(n)) '.mat'],'results_struct_ir', '-v7.3')
    disp('Done.')
end    
delete(wb);
%%
cw_vec_sp = results_struct_ir.sim_struct_sp(1).metric_array(:,cw_index);
spec_vec_sp = results_struct_ir.sim_struct_sp(1).metric_array(:,spec_index);
sh_vec_sp = results_struct_ir.sim_struct_sp(1).metric_array(:,sharpness_index);
sh_ft_sp = sh_vec_sp>=0;
sh_norm_vec_sp = results_struct_ir.sim_struct_sp(1).metric_array(:,sharpness_norm_index);
ir_vec_sp = results_struct_ir.sim_struct_sp(1).metric_array(:,ir_index);

%%
cw_vec_sh = results_struct_ir.sim_struct_sh(1).metric_array(:,cw_index);
sh_norm_vec_sh = results_struct_ir.sim_struct_sh(1).metric_array(:,sharpness_norm_index);
ir_vec_sh = results_struct_ir.sim_struct_sh(1).metric_array(:,ir_index);
%%
cw_vec_ir = results_struct_ir.sim_struct_ir(1).metric_array(:,cw_index);
phi_vec_ir = results_struct_ir.sim_struct_ir(1).metric_array(:,phi_index);
spec_vec_ir = results_struct_ir.sim_struct_ir(1).metric_array(:,spec_index);
spec_full_vec_ir = results_struct_ir.sim_struct_ir(1).metric_array(:,spec_full_index);
sh_vec_ir = results_struct_ir.sim_struct_ir(1).metric_array(:,sharpness_index);
sh_norm_vec_ir = results_struct_ir.sim_struct_ir(1).metric_array(:,sharpness_norm_index);
r_vec_ir = results_struct_ir.sim_struct_ir(1).metric_array(:,rate_index);
sh_vec_ir_norm = sh_vec_ir ./ (r_vec_ir.*(1-r_vec_ir));
sh_ft_ir = sh_vec_ir>=0;
ir_vec_ir = results_struct_ir.sim_struct_ir(1).metric_array(:,ir_index);

cw_bins = linspace(0,6,16);
cw_axis = cw_bins(1:end-1) + diff(cw_bins)/2;
max_ir_spec = NaN(1,length(cw_bins)-1);
max_ir_spec_full = NaN(1,length(cw_bins)-1);
max_ir_sharp = NaN(1,length(cw_bins)-1);
max_ir = NaN(1,length(cw_bins)-1);

for c = 1:length(cw_bins)-1
    cw_filter_ir = cw_vec_ir>=cw_bins(c) & cw_vec_ir<cw_bins(c+1);
%     cw_filter_sp = cw_vec_sp>=cw_bins(c) & cw_vec_sp<cw_bins(c+1);
%     cw_filter_sh = cw_vec_sh>=cw_bins(c) & cw_vec_sh<cw_bins(c+1);
    
    % get global max
    max_ir(c) = nanmax(ir_vec_ir(cw_filter_ir));
    ir_99 = prctile(ir_vec_ir(cw_filter_ir),99);
    
    % get specificity max
%     spec_val = prctile(spec_vec_sp(cw_filter_sp),98);
    max_ir_spec(c) = nanmean(spec_vec_ir(cw_filter_ir&ir_vec_ir>=ir_99));
    max_ir_spec_full(c) = nanmean(spec_full_vec_ir(cw_filter_ir&ir_vec_ir>=ir_99));
%     sharp_val = prctile(sh_norm_vec_sh(cw_filter_sh),98);
    max_ir_sharp(c) = nanmean(sh_norm_vec_ir(cw_filter_ir&ir_vec_ir>=ir_99));
end  
% %% Sharpness vs CW
% metric_indices_sh = [cw_bins sharpness_index];
% results_struct_sh = struct;
% for n = 1:length(nStateVec)-1
%     tic
%     [results_struct_sh(n).sim_info, results_struct_sh(n).sim_struct] = param_sweep_multi_v2(metric_indices_sh,functionPathCell{n}, sweep_options{:},...
%                                               'half_max_flag',true,'equilibrium_flag',true,'simType','General','nStates',nStateVec(n));
%     results_struct_sh(n).functionPath =  functionPathCell{n};
%     toc;    
% end    
% 
% % save output
% saveName = [metric_names{metric_indices_sh(1)} '_vs_' metric_names{metric_indices_sh(2)}];
% save([OutPath saveName '_eq.mat'],'results_struct_sh', '-v7.3')