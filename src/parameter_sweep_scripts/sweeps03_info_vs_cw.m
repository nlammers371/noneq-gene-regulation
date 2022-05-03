% script to call core parameter sweep function to examine tradeoffs between
% different transcritional input/output behaviors

clear 
close all
addpath(genpath('../utilities/'))

% as we add additional reaction, this will multiply the state space by
% a factor of 2, since each state from smaller "base" network can coincide
% with the new factor being (1) engaged or (2) disengaged

% set basic parameters
n_bs_vec = 1:5;
n_g_vec = 1:5;
ns_flag = 1;
ds_flag = 1;

% Set Dropbox directory
DropboxFolder = 'S:\Nick\Dropbox\Nonequilibrium\Nick\SweepOutput';
writePath = [DropboxFolder filesep 'sweeps03_info_vs_cw' filesep];
mkdir(writePath);

% this contains paths used to address correct functions
functionPathCell = cell(length(n_g_vec),length(n_bs_vec));
saveNameCell = cell(length(n_g_vec),length(n_bs_vec));

% get metric names
[~,~,metric_names] = calculateMetricsNumeric_v3([]);

% generate path to save metric functions 
subfolderName = 'numeric';
sourcePath = handlePathOptions(subfolderName);


for m = 1:length(n_g_vec)
    for n = 1:length(n_bs_vec) 

        % add path to functions        
        readName = ['s' sprintf('%02d',n_bs_vec(n)) '_ns00_g'...
                  sprintf('%02d',n_g_vec(m)) '_cw' num2str(ns_flag)];
        functionPath = [sourcePath readName filesep];
        functionPathCell{m,n} = functionPath;
        saveNameCell{m,n} = readName;

    end
end    
                      
% set sim options
sweep_options = {'n_sim',500,'n_seeds',15,'n_iters_max',50, 'numerical_precision',10, 'useParpool',1};
%%   
rate_index = find(strcmp(metric_names,'ProductionRate'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
precision_index = find(strcmp(metric_names,'Precision'));
phi_index = find(strcmp(metric_names,'Phi'));    
tau_index = find(strcmp(metric_names,'TauCycle'));  
cw_index = find(strcmp(metric_names,'CW'));  
ir_index = find(strcmp(metric_names,'IR'));
        
% paramBounds = repmat([-5; 5],1,12);
% paramBounds(1,6) = 0;
% paramBounds(2,7) = 0;
% results_struct = struct;
% iter = 1;
for equilibrium_flag = 0
    for m = 2:4
        if m == 2
            bs_list = 3:5;
        elseif m == 3
            bs_list = 3:4;
        elseif m == 4 
            bs_list = 3;
        end
        for n = bs_list
            % call sweep function
            tic
            [sim_info, sim_results, sim_results_short] = ...
                                param_sweep_multi_v3([cw_index ir_index],...
                                functionPathCell{m,n}, sweep_options{:},...
                                'equilibrium_flag',equilibrium_flag,'TauCycleTime',1,...
                                'downsample_output',ds_flag);

            if ds_flag
                sim_results = sim_results_short;
            end            
            suffix = '';
            if equilibrium_flag
                suffix = '_eq';
            end
            % save
            disp('saving...')
            saveName = saveNameCell{m,n};
            save([writePath 'sweep_info_' saveName suffix '.mat'],'sim_info')
            save([writePath 'sweep_results_' saveName suffix '.mat'],'sim_results', '-v7.3')

            toc;    
    %         iter = iter + 1;
        end
    end    
end
%%
[ir_max, ir_i] = nanmax(sim_results.metric_array(:,ir_index))
rates_opt = log10(sim_results.rate_array(ir_i,:))

% close all
% 
% figure;
% hold on
% for n = 1
%     
%     phi_index = find(strcmp(metric_names,'Phi'));
%     ir_index = find(strcmp(metric_names,'DecisionRateNorm'));
%     sharpness_index = find(strcmp(metric_names,'Sharpness'));
%     
%     % concatenate arrays
% %     metric_array = [];
% %     rate_array_long = [];
% %     for ns_flag = 1:length(results_struct(n).sim_struct)
% %         metric_array = vertcat(metric_array,results_struct(n).sim_struct(ns_flag).metric_array);
% %         rate_array_long = vertcat(rate_array_long,results_struct(n).sim_struct(ns_flag).rate_array);
% %     end
% %     results_struct(n).sim_struct.metric_array = metric_array;
% %     results_struct(n).sim_struct.rate_array = rate_array_long;
%     phi_filter = metric_array(:,phi_index)>=-20;
%     tau_filter = metric_array(:,tau_index)<=10*3600 & metric_array(:,tau_index)>=0;
% %     scatter(exp(metric_array(:,phi_index)),metric_array(:,ir_index))    
%     scatter(exp(metric_array(phi_filter&tau_filter,phi_index)),metric_array(phi_filter&tau_filter,ir_index))    
% end
% xlim([1e-4 1e4])
% % ylim([0 0.03])
% set(gca,'xscale','log');
% 
% %%
% n = 1;
% w_flags = contains(results_struct(n).sim_info.sweepVarStrings,'w');
% sweepFlags = results_struct(n).sim_info.sweepFlags;
% 
% tau_index = find(strcmp(metric_names,'CycleTime'));
% phi_multiplier = (1*results_struct(n).sim_struct.metric_array(:,phi_index)>=0);
% tau_vec = results_struct(n).sim_struct.metric_array(:,tau_index);
% tau_multiplier = (1*(tau_vec>=0&tau_vec<=3600*10));
% [ir_max,ir_i] = nanmax(results_struct(n).sim_struct.metric_array(:,ir_index))
% tau_c = results_struct(n).sim_struct.metric_array(ir_i,tau_index);
% rates_raw = results_struct(n).sim_struct.rate_array(ir_i,:)
% rates_adjusted = log10(rates_raw);
% % rates_adjusted(sweepFlags&~w_flags) = rates_adjusted(sweepFlags&~w_flags).*tau_c/300;
% % rates_log = log10(rates_adjusted)
% 
% %%
% [~,~,metric_names] = calculateMetricsNumeric_v3([]);
% phi_index = find(strcmp(metric_names,'Phi'));
% ir_index = find(strcmp(metric_names,'DecisionRateNorm'));
% rate_index = find(strcmp(metric_names,'Production Rate'));
% 
% 
% phi_vec = results_struct(2).metric_array(:,phi_index);
% ir_vec = results_struct(2).metric_array(:,ir_index);
% r_vec = results_struct(2).metric_array(:,rate_index);
% tau_vec =results_struct(2).metric_array(:,tau_index);
% 
% [ir_max,ir_i] = nanmax(ir_vec.*(1*phi_vec<=5e2).*(1*phi_vec>=0).*(1*tau_vec<=300))
% 
% best_parameters = results_struct(2).rate_array_long(ir_i,:)
% 
% %% Look at architectures of sharp networks
% state_prob_cell = cell(1,length(results_struct));
% rate_cell = cell(1,length(results_struct));
% opt_index_vec = NaN(1,length(results_struct));
% opt_sharpness_vec = NaN(1,length(results_struct));
% 
% for n = 1:length(results_struct)
%     sharp_filter = results_struct(n).sim_struct(1).metric_array(:,sharpness_index)>=0;
%     % get index of sharpest network
%     [opt_sharpness_vec(n) ,opt_index_vec(n)] = nanmax(results_struct(n).sim_struct(1).metric_array(:,sharpness_index).*sharp_filter);
%     % generate network
%     opt_params = results_struct(n).sim_struct(1).rate_array(opt_index_vec(n),:);
%     valCell = mat2cell(opt_params,size(opt_params,1),ones(1,size(opt_params,2)));   
%     % map to appropriate functions
%     rmpath(genpath('../utilities/metricFunctions/'));
%     addpath(genpath(results_struct(n).sim_info.functionPath));
%     % get rate array
%     rate_cell{n} = RSymFun(valCell{:});
%     % get probs
%     state_prob_cell{n} = calculate_ss_num(rate_cell{n},10);      
% end    
% 
