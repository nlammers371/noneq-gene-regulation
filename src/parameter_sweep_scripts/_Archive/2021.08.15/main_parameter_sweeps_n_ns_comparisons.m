% script to call core parameter sweep function to examine tradeoffs between
% different network behaviors

clear 
close all
addpath(genpath('../utilities/'))

[~,~,metric_names] = calculateMetricsNumeric_v3([]);
nStateVec = [4 8];



for n = 1:length(nStateVec)
    if n > 1
      f_type = 'numeric';
      
      % generate path to save metric functions 
      functionPath = handlePathOptions(subfolderName);
      
      functionPathCell(n) = {['../utilities/metricFunctions/n' num2str(nStateVec(n)) '_OR_NS_NUM/']};        
    else
      f_type = 'symbolic';
      functionPathCell(n) = {['../utilities/metricFunctions/n4_OR/']};
    end
end
                         
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
sweep_options = {'n_sim',1,'n_seeds',5,'n_iters_max',50,'numCalcFlag',1,'numerical_precision',10};
%%

% results_struct = struct;
for n = 1:length(nStateVec)
    tic
    [results_struct(n).sim_info, results_struct(n).sim_struct] = param_sweep_multi_v2([phi_index ir_index],functionPathCell{n}, sweep_options{:},...
                                              'half_max_flag',false,'equilibrium_flag',false,'simType','General','nStates',nStateVec(n));
    toc;    
%     tic
%     [results_struct(n).sim_info_eq, results_struct(n).sim_struct_eq] = param_sweep_multi_v2([rate_index sharpness_index],functionPathCell{n}, sweep_options{:},...
%                                               'half_max_flag',false,'equilibrium_flag',true,'simType','General','nStates',nStateVec(n));
%     toc;
end    

%%
% close all

figure;
hold on
for n = 1:2%3:-1:1%fliplr(1:4)
    sharp_filter = results_struct(n).sim_struct(1).metric_array(:,sharpness_index)>=0;
    scatter(results_struct(n).sim_struct(1).metric_array(sharp_filter,phi_index),...
      results_struct(n).sim_struct(1).metric_array(sharp_filter,ir_index))
    
end

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

