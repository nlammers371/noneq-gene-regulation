% script to call core parameter sweep function to examine tradeoffs between
% different transcritional input/output behaviors

% clear 
% close all
addpath(genpath('../utilities/'))

% as we add additional reaction, this will multiply the state space by
% a factor of 2, since each state from smaller "base" network can coincide
% with the new factor being (1) engaged or (2) disengaged

% set basic parameters
n_bs_vec = 1:2;
n_g_vec = 1:2;
ns_flag_vec = 1;
% n_bs_cw = [6 12 20 30 42]/2;

% nStateArray = [6 18; 12 36];
% subfolderName = 'numeric\unique';
% nStateArray = [6 12; 12 24];
% subfolderName = 'numeric\unique_g';
% nStateArray = [6 18; 12 36];
% subfolderName = 'numeric\unique_a';
subfolderNameCell = {'numeric\unique','numeric\unique_a','numeric\unique_g','numeric'};

functionPathCell = cell(length(n_g_vec),length(n_bs_vec),length(subfolderNameCell));
% iter = 1;

% get metric names
[~,~,metric_names] = calculateMetricsNumeric_v3([]);

% generate path to save metric functions 
ns_flag = 1;


for i = 1:length(subfolderNameCell)
    writePath = handlePathOptions(subfolderNameCell{i});
    for m = 1:length(n_g_vec)
        for n = 1:length(n_bs_vec) 
            
            % add path to functions
%             nStates = nStateArray(m,n,i);
            functionPath = [writePath 's' sprintf('%02d',n_bs_vec(n)) '_ns00_g' sprintf('%02d',n_g_vec(m)) '_cw' num2str(ns_flag) filesep];
            functionPathCell{m,n,i} = functionPath;
    
        end
    end    
end
                         
% set sim options
sweep_options = {'n_sim',24,'n_seeds',10,'n_iters_max',50,'numerical_precision',10,'useParpool',1};
%%
rate_index = find(strcmp(metric_names,'Production Rate'));    
sharpness_index = find(strcmp(metric_names,'Sharpness'));
cw_index = find(strcmp(metric_names,'CW'));
phi_index = find(strcmp(metric_names,'Phi'));    
ir_index = find(strcmp(metric_names,'DecisionRateNorm'));
tau_index = find(strcmp(metric_names,'CycleTime'));

% results_struct = struct;
% iter = 1;
for i =1% CW or not?
    for m = 2% % how many general reactions?
        for n = 1:2 % how many binding sites?

            % add path to functions
%             nStates = nStateArray(m,n);

            % call sweep function
            tic
            [results_struct(iter).sim_info, results_struct(iter).sim_struct] = ...
                                param_sweep_multi_v3([phi_index ir_index],...
                                functionPathCell{m,n,i}, sweep_options{:},...
                                'half_max_flag',false,'cw',1,'equilibrium_flag',false,'TauCycleTime',3600);
            toc;    
            iter = iter + 1;
        end
    end    
end
%%
% close all
% 
figure;
hold on
for n = [14 10]%-1:1%6:-1:1%length(nStateVec):-1:1\
    
    % concatenate arrays
    metric_array_long = [];
    rate_array_long = [];
    for ns_flag = 1:length(results_struct(n).sim_struct)
        metric_array_long = vertcat(metric_array_long,results_struct(n).sim_struct(ns_flag).metric_array);
        rate_array_long = vertcat(rate_array_long,results_struct(n).sim_struct(ns_flag).rate_array);
    end
    results_struct(n).metric_array_long = metric_array_long;
    results_struct(n).rate_array_long = rate_array_long;
    
    metric_array_long = results_struct(n).metric_array_long;
    
    phi_filter = metric_array_long(:,phi_index)>=-20;
    tau_filter =  metric_array_long(:,tau_index)<=10*3600 & metric_array_long(:,tau_index)>=0;
%     scatter(exp(metric_array_long(:,phi_index)),metric_array_long(:,ir_index))    
    scatter(exp(metric_array_long(phi_filter&tau_filter,phi_index)),metric_array_long(phi_filter&tau_filter,ir_index))    
end
xlim([1e-5 1e4])
ylim([1e-8 2e1])
set(gca,'yscale','log');
set(gca,'xscale','log');

%%
n = 7;
w_flags = contains(results_struct(n).sim_info.sweepVarStrings,'w');
sweepFlags = results_struct(n).sim_info.sweepFlags;

tau_index = find(strcmp(metric_names,'CycleTime'));
cw_filter = 1*(results_struct(n).metric_array_long(:,rate_index)<=0.51&results_struct(n).metric_array_long(:,rate_index)>=0.49);
phi_multiplier = (1*results_struct(n).metric_array_long(:,phi_index)>=-20);
tau_vec = results_struct(n).metric_array_long(:,tau_index);
tau_multiplier = (1*(tau_vec>=0&tau_vec<=3600*10));

[ir_max,ir_i] = nanmax(results_struct(n).metric_array_long(:,ir_index).*tau_multiplier.*phi_multiplier)

tau_c = results_struct(n).metric_array_long(ir_i,tau_index);
rates_raw = results_struct(n).rate_array_long(ir_i,:);
rates_adjusted = rates_raw;
rates_adjusted(sweepFlags&~w_flags) = rates_adjusted(sweepFlags&~w_flags).*tau_c/300;
% rates_log = log10(rates_adjusted)

%%
[~,~,metric_names] = calculateMetricsNumeric_v3([]);
phi_index = find(strcmp(metric_names,'Phi'));
ir_index = find(strcmp(metric_names,'DecisionRateNorm'));
rate_index = find(strcmp(metric_names,'Production Rate'));


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

