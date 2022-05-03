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
nStateVec = [4 8 16 32];
% numerical_precision = 5;
for n = 1:length(nStateVec)
    functionPathCell(n) = {['../utilities/metricFunctions/n' num2str(nStateVec(n)) '_OR_SPEC_NUM/']};    
end

OutPath = [DropboxPath 'parameter_sweeps_multi_right_site/'];
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
tau_index = find(strcmp(metric_names,'CycleTime'));

% set sim options
sweep_options = {'n_sim',1,'n_seeds',10,'n_iters_max',50,'numCalcFlag',1,'numerical_precision',10,'useParpool',0,'NumWorkers',4};

%% IR vs CW for each motif
% metric_indices_ir = [cw_index ir_index];

totalRuns = length(nStateVec)*3;
iter = 1;
sweep_results = struct;
for n = 1%:length(nStateVec)
    
    tic
    %%%%%%%%%%%%%%%
    % transcription rate vs sharpness 
    % NEQ
    [sweep_results(iter).sim_info_sh, sweep_results(iter).sim_struct_sh] = param_sweep_multi_v3([sharpness_index precision_index],...
                                              functionPathCell{n}, sweep_options{:},'cycleTime',100,...
                                              'half_max_flag', false,'equilibrium_flag',false,'nStates',nStateVec(n));
    sweep_results(iter).sh_opt_indices = [cw_index ir_index];
    
    sweep_results(iter).functionPath =  functionPathCell{n};
   
    % EQ
%     [sweep_results(iter).sim_info_sh_eq, sweep_results(iter).sim_struct_sh_eq] = param_sweep_multi_v3([sharpness_index precision_index],...
%                                               functionPathCell{n}, sweep_options{:},'cycleTime',100,...
%                                               'half_max_flag', false,'equilibrium_flag',true,'nStates',nStateVec(n));
   
    toc
    iter = iter + 1;
end    
disp(['Saving results...'])
save([OutPath 'sweep_results.mat'],'sweep_results', '-v7.3')
disp('Done.')
%%
% first let's remove networks that are that violate parameter bounds once
% cycle time normalization is applied
% tau_target = 100;
% tau_vec = sweep_results(1).sim_struct_sh(1).metric_array(:,tau_index);
% rates_adjusted = sweep_results(1).sim_struct_sh(1).rate_array ./ tau_vec .*tau_target;
% 
% paramBounds = sweep_results(1).sim_info_sh(1).paramBounds;
% sweepFlags = sweep_results(1).sim_info_sh(1).sweepFlags;
% 
% oob_flags = max(log10(rates_adjusted(:,sweepFlags))<paramBounds(1,sweepFlags)...
%                 |log10(rates_adjusted(:,sweepFlags))>paramBounds(2,sweepFlags),[],2);
%               
% %%              
%               

close all
figure;
hold on
scatter(sweep_results(1).sim_struct_sh(1).metric_array(:,sharpness_index),...
        16*100*exp(sweep_results(1).sim_struct_sh(1).metric_array(:,precision_index)).^2)
      
% scatter(sweep_results(1).sim_struct_sh_eq(1).metric_array(:,sharpness_index),...
%         16*100*exp(sweep_results(1).sim_struct_sh_eq(1).metric_array(:,precision_index)).^2)      

%%
[max_sharp, mi_index] = nanmax(sweep_results(1).sim_struct_sh(1).metric_array(:,sharpness_index));
max_rates = sweep_results(1).sim_struct_sh(1).rate_array(mi_index,:);
max_binding_rates = max_rates(sweep_results(1).sim_info_sh(1).bindingFlags==1)
max_unbinding_rates = max_rates(sweep_results(1).sim_info_sh(1).unbindingFlags==1)

% %%
% % close all
% % 
% % n_states_to_plot = [8];
% sharpness_vs_rate = figure;
% hold on
% for n  = 4:-1:1%length(sweep_results)
%     rate_vec = [];
%     sharpness_vec = [];
%     for s = 1:length(sweep_results(n).sim_struct_sh)
%         sharpness_vec = vertcat(sweep_results(n).sim_struct_sh(s).metric_array(:,sharpness_index),sharpness_vec);
%         rate_vec = vertcat(sweep_results(n).sim_struct_sh(s).metric_array(:,rate_index),rate_vec);
%     end
% %     plot_filter = sharpness_vec>=0;
%     scatter(rate_vec,sharpness_vec)
% end
%   
%   
%   
%   
