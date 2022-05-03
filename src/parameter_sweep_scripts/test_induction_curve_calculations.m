% script to call core parameter sweep function to examine tradeoffs between
% different network behaviors
clear 
close all
addpath(genpath('../utilities/'))

[~,~,metric_names] = calculateMetricsNumeric_v3([]);
nStates = 18;
folder_prefix = 's02_ns00_g01_cw0';
% functionPathCell = cell(size(nStates));
% P:\Nick\projects\noneq-transcription\src\utilities\metricFunctions\numeric\n018_s02_ns00_g02_cw1
f_type = 'numeric'; 
% generate path to save metric functions 
topFolderPath = handlePathOptions(f_type);
% generate subfolder name
subFolderName = [folder_prefix filesep];
functionPath = [topFolderPath subFolderName];                            
                         
% set sim options
sweep_options = {'n_sim',0,'n_seeds',10,'n_iters_max',50,'numerical_precision',10};

% run dummy trial to get sim info vec

sharp_index = find(strcmp(metric_names,'Sharpness'));    
prec_index = find(strcmp(metric_names,'Precision'));

% call sweep function
tic
[simInfo, sim_struct] = ...
                    param_sweep_multi_v3([sharp_index prec_index],...
                    functionPath, sweep_options{:},...
                    'half_max_flag',false,'equilibrium_flag',false);
                  
% set param values for static variables
n_samples = 500;
sweepFlags = simInfo.sweepFlags;
param_array = NaN(n_samples,size(sweepFlags,2));
defaultValues = simInfo.defaultValues;
paramBounds = simInfo.paramBounds;
simInfo.rep_size = 10;%*simInfo.min_points_per_bin
% set values for static 
param_array(:,~sweepFlags) = repmat(defaultValues(~sweepFlags),size(param_array,1),1);

% generate initial array of samples
lb_array = repmat(paramBounds(1,sweepFlags)/4,n_samples,1);
ub_array = repmat(paramBounds(2,sweepFlags)/4,n_samples,1);

param_array(:,sweepFlags) = sample_rates_multi_v2(lb_array,ub_array,[],4,simInfo); 
param_array(:,2) = 1;
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));
%% calculate stuff
c_vec = logspace(-3,3,150)';

pd_rate_full = NaN(n_samples,1);
pd_rate = NaN(n_samples,1);

var_full = NaN(n_samples,1);
var_short = NaN(n_samples,1);

numerical_precision = 2;

for n = 1%:n_samples
    param_temp = repmat(param_array(n,:),length(c_vec),1);
    param_temp(:,1) = c_vec;         

    % get full rate matrix   
%     Q_num_full = RSymFunFull(valCellCS{:});
    pd_rate_temp = NaN(size(c_vec));
    var_temp = NaN(size(c_vec));
    tic
    for c = 1:length(c_vec)
        valCellCS = mat2cell(param_temp(c,:),size(param_temp(c,:),1),ones(1,size(param_temp(c,:),2)));   
        Q_num = RSymFun(valCellCS{:});
    
    % state probabilities
%     ss_full = calculate_ss_num(Q_num_full,numerical_precision);  
%     pd_rate_full(n) = sum(ss_full(simInfo.activeStateFilterFull));
     
        ss_short = calculate_ss_num(Q_num,numerical_precision);  
        pd_rate_temp(c) = sum(ss_short(simInfo.activeStateFilter));
        Z_num_short = calculate_Z_matrix(Q_num,ss_short,numerical_precision);
        var_temp(c) = calculate_var_num(Q_num,ss_short,simInfo.activeStateFilter,Z_num_short,numerical_precision);
    end
    toc
%     % variance
%     Z_num_full = calculate_Z_matrix(Q_num_full,ss_full,numerical_precision);
%     var_full(n) = calculate_var_num(Q_num_full,ss_full,simInfo.activeStateFilterFull,Z_num_full,numerical_precision);
%     
    tic
    Z_num_short = calculate_Z_matrix(Q_num,ss_short,numerical_precision);
    toc
    tic
    [TauOn,TauOff,TauCycle] = calculate_tau_num(Q_num,ss_short,simInfo.activeStateFilter,numerical_precision);
    toc
%     var_short(n) = calculate_var_num(Q_num,ss_short,simInfo.activeStateFilter,Z_num_short,numerical_precision);
    
end    

%%
close all
figure;
scatter(pd_rate,pd_rate_full)

figure;
scatter(var_short,var_full)
set(gca,'yscale','log')
set(gca,'xscale','log')