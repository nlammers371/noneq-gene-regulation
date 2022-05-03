% script to call core parameter sweep function to examine tradeoffs between
% different network behaviors
clear 
close all
addpath(genpath('../utilities/'))

% [~,~,metric_names] = calculateMetricsNumeric_v3([]);
nStateVec = [4 6 8 10 12];
folder_prefix = '_ns00_g01';
functionPathCell = cell(size(nStateVec));

for n = 1:length(nStateVec)
    if true%n > 1
      f_type = 'numeric'; 
    else
      f_type = 'symbolic';
    end
      % generate path to save metric functions 
      topFolderPath = handlePathOptions(f_type);
      % generate subfolder name
      subFolderName = ['n' sprintf('%03d',nStateVec(n)) '_s' sprintf('%02d',n) folder_prefix filesep];
      functionPathCell(n) = {[topFolderPath subFolderName]};                            
end
                         
% set sim options
sweep_options = {'n_sim',10,'n_seeds',10,'n_iters_max',50,'numerical_precision',10};
%%

results_struct = struct;
for n = 3%:length(nStateVec)
    % get metric indices
    if false%n == 1
      [~,~,metric_names] = calculateMetricsSym_v2([]);        
    else
      [~,~,metric_names] = calculateMetricsNumeric_v3([]);
    end
    
    phi_index = find(strcmp(metric_names,'Phi'));    
    ir_index = find(strcmp(metric_names,'DecisionRateNorm'));

    % call sweep function
    tic
    [results_struct(n).sim_info, results_struct(n).sim_struct] = ...
                        param_sweep_multi_v3([phi_index ir_index],...
                        functionPathCell{n}, sweep_options{:},...
                        'half_max_flag',false,'equilibrium_flag',false,'nStates',nStateVec(n));
    toc;    

end    

%%
close all

figure;
hold on
for n = 3%length(nStateVec):-1:1
    if false
        [~,~,metric_names] = calculateMetricsSym_v2([]);        
    else
        [~,~,metric_names] = calculateMetricsNumeric_v3([]);
    end
    phi_index = find(strcmp(metric_names,'Phi'));
    ir_index = find(strcmp(metric_names,'DecisionRateNorm'));
%     tau_index = find(strcmp(metric_names,'CycleTime'));
    sharpness_index = find(strcmp(metric_names,'Sharpness'));
    
    % concatenate arrays
    metric_array_long = [];
    rate_array_long = [];
    for i = 1:length(results_struct(n).sim_struct)
        metric_array_long = vertcat(metric_array_long,results_struct(n).sim_struct(i).metric_array);
        rate_array_long = vertcat(rate_array_long,results_struct(n).sim_struct(i).rate_array);
    end
    results_struct(n).metric_array_long = metric_array_long;
    results_struct(n).rate_array_long = rate_array_long;
    
    scatter(metric_array_long(:,phi_index),metric_array_long(:,ir_index))    
end
xlim([1e-4 1e4])
% ylim([0 0.03])
set(gca,'xscale','log');

%%
[~,~,metric_names] = calculateMetricsNumeric_v3([]);
phi_index = find(strcmp(metric_names,'Phi'));
ir_index = find(strcmp(metric_names,'DecisionRateNorm'));
rate_index = find(strcmp(metric_names,'Production Rate'));
tau_index = find(strcmp(metric_names,'CycleTime'));

phi_vec = results_struct(2).metric_array_long(:,phi_index);
ir_vec = results_struct(2).metric_array_long(:,ir_index);
r_vec = results_struct(2).metric_array_long(:,rate_index);
tau_vec =results_struct(2).metric_array_long(:,tau_index);

[ir_max,ir_i] = nanmax(ir_vec.*(1*phi_vec<=5e2).*(1*phi_vec>=0).*(1*tau_vec<=300))

best_parameters = results_struct(2).rate_array_long(ir_i,:)

%% Look at architectures of sharp networks
state_prob_cell = cell(1,length(results_struct));
rate_cell = cell(1,length(results_struct));
opt_index_vec = NaN(1,length(results_struct));
opt_sharpness_vec = NaN(1,length(results_struct));

for n = 1:length(results_struct)
    sharp_filter = results_struct(n).sim_struct(1).metric_array(:,sharpness_index)>=0;
    % get index of sharpest network
    [opt_sharpness_vec(n) ,opt_index_vec(n)] = nanmax(results_struct(n).sim_struct(1).metric_array(:,sharpness_index).*sharp_filter);
    % generate network
    opt_params = results_struct(n).sim_struct(1).rate_array(opt_index_vec(n),:);
    valCell = mat2cell(opt_params,size(opt_params,1),ones(1,size(opt_params,2)));   
    % map to appropriate functions
    rmpath(genpath('../utilities/metricFunctions/'));
    addpath(genpath(results_struct(n).sim_info.functionPath));
    % get rate array
    rate_cell{n} = RSymFun(valCell{:});
    % get probs
    state_prob_cell{n} = calculate_ss_num(rate_cell{n},10);      
end    

