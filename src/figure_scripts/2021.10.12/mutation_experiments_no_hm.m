% script to track information rate as a function of cw

clear 
close all
addpath(genpath('../utilities/'))

% set basic parameters
nStates = 6;
paramBounds = repmat([-10 ; 6],1,11); % constrain transition rate magnitude
[~,~,metric_names] = calculateMetricsSym_v2([]);

% specify function path
functionPath = ['../utilities/metricFunctions/n' num2str(nStates) '_OR/'];

% make sure we're linked to the appropriate function subfolder% make sure we're linked to the appropriate function subfolder
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));

% define save path

DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\manuscript\';
% DropboxFolder = 'S:\Nick\Dropbox\Nonequilibrium\Nick\manuscript\';

FigPath = [DropboxFolder 'observed_effects' filesep];
mkdir(FigPath);         

% get index of useful metrics
spec_index = find(strcmp(metric_names,'Specificity'));
sharp_right_norm_index = find(strcmp(metric_names,'SharpnessRightNorm'));
ir_index = find(strcmp(metric_names,'DecisionRateNorm'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
sharpness_norm_index = find(strcmp(metric_names,'SharpnessNormed'));
precision_index = find(strcmp(metric_names,'Precision'));
rate_index = find(strcmp(metric_names,'Production Rate'));
precision_right_index = find(strcmp(metric_names,'PrecisionRight'));
cw_index = find(strcmp(metric_names,'CW'));

% set sim options
sweep_options = {'n_seeds',5,'n_iters_max',50,'n_sim',10,'nStates',nStates};

% calculate sensitivity bound
alpha_factor = 100;

% specify plot options 
markerAlpha = 0.5; % marker transparency
markerSize = 75; % marker size
n_plot = 3e3;
bit_factor = log2(exp(1));

%% Run initial sweeps
cw_vec = 0:0.5:4;

% cw = 1e5; %
master_struct = struct;
for c = 1:length(cw_vec)
    tic
    [sim_info_neq, sim_struct_neq] = param_sweep_multi_v3([ir_index sharpness_index],functionPath,sweep_options{:},...
                                              'half_max_flag',false,'cw',10.^cw_vec(c),...
                                              'equilibrium_flag',false,'specFactor',alpha_factor,'paramBounds',paramBounds);
                                            
    [sim_info_eq, sim_struct_eq] = param_sweep_multi_v3([ir_index sharpness_index],functionPath,sweep_options{:},...
                                              'half_max_flag',false,'cw',10.^cw_vec(c),...
                                              'equilibrium_flag',true,'specFactor',alpha_factor,'paramBounds',paramBounds);  
                                            
    toc
                                            
    % extract key pieces of info
    
    % eq
    metric_array_eq = vertcat(sim_struct_eq.metric_array);
    rate_array_eq = vertcat(sim_struct_eq.rate_array);
    
    sharp_vec_eq = metric_array_eq(:,sharpness_index);
    ir_vec_eq = metric_array_eq(:,ir_index);
    ir_99_eq = prctile(ir_vec_eq,99.5);
    sharp_99_eq = prctile(sharp_vec_eq,99);
    optimal_ids_eq = find(ir_vec_eq>=ir_99_eq & sharp_vec_eq>=0);%>=sharp_99_eq);

    master_struct(c).metric_array_eq = metric_array_eq;
    master_struct(c).rate_array_eq = rate_array_eq;
    master_struct(c).opt_rates_eq = rate_array_eq(optimal_ids_eq,:);
    master_struct(c).opt_ir_eq = ir_vec_eq(optimal_ids_eq);
    master_struct(c).opt_s_eq = sharp_vec_eq(optimal_ids_eq);
    master_struct(c).opt_ids_eq = optimal_ids_eq;
    
    master_struct(c).cw_val_eq = repmat(cw_vec(c),length(optimal_ids_eq),1);
    
    % neq (global)
    metric_array_neq = vertcat(sim_struct_neq.metric_array);
    rate_array_neq = vertcat(sim_struct_neq.rate_array);
    
    
    sharp_vec_neq = metric_array_neq(:,sharpness_index);
    ir_vec_neq = metric_array_neq(:,ir_index);
    ir_99_neq = prctile(ir_vec_neq,99.5);
    sharp_99_neq = prctile(sharp_vec_neq,99);
    optimal_ids_neq = find(ir_vec_neq>=ir_99_neq & sharp_vec_neq>=0);

    master_struct(c).metric_array_neq = metric_array_neq;
    master_struct(c).rate_array_neq = rate_array_neq;
    master_struct(c).opt_rates_neq = rate_array_neq(optimal_ids_neq,:);
    master_struct(c).opt_ir_neq = ir_vec_neq(optimal_ids_neq);
    master_struct(c).opt_s_neq = sharp_vec_neq(optimal_ids_neq);
    master_struct(c).opt_ids_neq = optimal_ids_neq;
                    
    % store cw values
    master_struct(c).cw_val_neq = repmat(cw_vec(c),length(optimal_ids_neq),1);
end
% %%
% tic
% [sim_info_neq, sim_struct_neq1] = param_sweep_multi_v3([ir_index sharpness_index],functionPath,sweep_options{:},...
%                                               'half_max_flag',false,'cw',10.^cw_vec(end),...
%                                               'equilibrium_flag',false,'specFactor',alpha_factor,'n_seeds',10,'paramBounds',paramBounds);
%                                             
% [sim_info_neq, sim_struct_neq2] = param_sweep_multi_v3([ir_index sharpness_index],functionPath,sweep_options{:},...
%                                               'half_max_flag',true,'cw',10.^cw_vec(end),...
%                                               'equilibrium_flag',false,'specFactor',alpha_factor,'n_seeds',10,'paramBounds',paramBounds);                                            
% 
% [sim_info_eq, sim_struct_eq1] = param_sweep_multi_v3([ir_index sharpness_index],functionPath,sweep_options{:},...
%                                               'half_max_flag',true,'cw',10.^cw_vec(end),...
%                                               'equilibrium_flag',true,'specFactor',alpha_factor,'n_seeds',10,'paramBounds',paramBounds);
%                                             
% [sim_info_eq, sim_struct_eq2] = param_sweep_multi_v3([ir_index sharpness_index],functionPath,sweep_options{:},...
%                                               'half_max_flag',false,'cw',10.^cw_vec(end),...
%                                               'equilibrium_flag',true,'specFactor',alpha_factor,'n_seeds',10,'paramBounds',paramBounds);                                            
% toc
% 
% %%
% test_eq1 = vertcat(sim_struct_eq1.metric_array);
% test_eq2 = vertcat(sim_struct_eq2.metric_array);
% 
% test_neq1 = vertcat(sim_struct_neq1.metric_array);
% test_neq2 = vertcat(sim_struct_neq2.metric_array);
% 
% ir_max_eq1 = nanmax(test_eq1(:,ir_index))
% ir_max_eq2 = nanmax(test_eq2(:,ir_index))
% ir_max_neq1 = nanmax(test_neq1(:,ir_index))
% ir_max_neq2 = nanmax(test_neq2(:,ir_index))
%%
close all
figure;
scatter(vertcat(master_struct.cw_val_eq),vertcat(master_struct.opt_ir_eq))
hold on
scatter(vertcat(master_struct.cw_val_neq),vertcat(master_struct.opt_ir_neq))

figure;
scatter(vertcat(master_struct.cw_val_eq),vertcat(master_struct.opt_s_eq))
hold on
scatter(vertcat(master_struct.cw_val_neq),vertcat(master_struct.opt_s_neq))

%% Simulate perturbations to binding site of differing levels of severity
% add functions to path
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));

% set basic parameters
n_perturbations = 1e2;
m_factor = logspace(0,log10(alpha_factor),n_perturbations);
alpha_adjusted = alpha_factor ./ m_factor;
% c_vec = logspace(-3,6,1e4);
% c_index = find(c_vec==1);


sharp_rate_array_eq = repmat(vertcat(master_struct.opt_rates_eq), 1,1,length(m_factor));
sharp_rate_array_eq(:,sim_info_eq.unbindingFlags==1,:) = sharp_rate_array_eq(:,sim_info_eq.unbindingFlags==1,:) .* reshape(m_factor,1,1,[]);
sharp_rate_array_eq(:,sim_info_eq.b_index,:) = repmat(reshape(alpha_adjusted,1,1,[]),size(sharp_rate_array_eq,1),1,1);

sharp_rate_array_neq = repmat(vertcat(master_struct.opt_rates_neq),1,1, length(m_factor));
sharp_rate_array_neq(:,sim_info_neq.unbindingFlags==1,:) = sharp_rate_array_neq(:,sim_info_eq.unbindingFlags==1,:) .* reshape(m_factor,1,1,[]);
sharp_rate_array_neq(:,sim_info_neq.b_index,:) = repmat(reshape(alpha_adjusted,1,1,[]),size(sharp_rate_array_neq,1),1,1);


%%
% calculate half-max point and sharpness 
cw_vec_eq = vertcat(master_struct.cw_val_eq);
cw_vec_neq = vertcat(master_struct.cw_val_neq);

sharpness_vec_eq = NaN(size(sharp_rate_array_eq,1),size(sharp_rate_array_eq,3));
sharpness_vec_neq = NaN(size(sharp_rate_array_neq,1),size(sharp_rate_array_neq,3));

rate_vec_eq = NaN(size(sharp_rate_array_eq,1),size(sharp_rate_array_eq,3));
rate_vec_neq = NaN(size(sharp_rate_array_neq,1),size(sharp_rate_array_neq,3));

for p = 1:size(sharp_rate_array_eq,3)
  
    % extract rates
    eq_rates = sharp_rate_array_eq(:,:,p);    
    neq_rates = sharp_rate_array_neq(:,:,p);        
    
    % convert to cell array and run calculations
    paramCellEq = mat2cell(eq_rates,size(eq_rates,1),ones(1,size(eq_rates,2)));
    paramCellNeq = mat2cell(neq_rates,size(neq_rates,1),ones(1,size(neq_rates,2)));
    
    % production rate 
    ProductionRateEq = productionRateFunction(paramCellEq{:});
    ProductionRateNeq = productionRateFunction(paramCellNeq{:});
    
    % sharpness
    SharpnessEq = sharpnessFunction(paramCellEq{:});
    SharpnessNeq = sharpnessFunction(paramCellNeq{:});    
    
    % record 
    sharpness_vec_eq(:,p) = SharpnessEq;
    sharpness_vec_neq(:,p) = SharpnessNeq;
    
    rate_vec_eq(:,p) = ProductionRateEq;
    rate_vec_neq(:,p) = ProductionRateNeq;
    
end

% %% Let's get a broader view of what is possible at equilibrium
% [sim_info_eq, sim_struct_eq] = param_sweep_multi_v3([cw_index rate_index],functionPath,sweep_options{:},...
%                                               'half_max_flag',false,...
%                                               'equilibrium_flag',true,'specFactor',alpha_factor); 
% 
% %%                                            
% rates_eq = vertcat(sim_struct_eq.rate_array);
% metrics_eq = vertcat(sim_struct_eq.metric_array);
% cw_vec_eq = metrics_eq(:,cw_index);
% 
% m_factor_index = 51; 
% rates_eq_mutated = rates_eq;
% rates_eq_mutated(:,sim_info_eq.unbindingFlags==1) = rates_eq_mutated(:,sim_info_eq.unbindingFlags==1) * m_factor(m_factor_index);
% rates_eq_mutated(sim_info_eq.b_index,:) = alpha_adjusted(m_factor_index);
% 
% paramCellEqMut = mat2cell(rates_eq_mutated,size(rates_eq_mutated,1),ones(1,size(rates_eq_mutated,2)));
% sharp_vec_eq_orig = metrics_eq(:,sharpness_index);
% sharp_vec_eq_mut = sharpnessFunction(paramCellEqMut{:});
% pos_filter = sharp_vec_eq_mut>=0 & sharp_vec_eq_orig >= 0;
% 
% fold_change_vec = sharp_vec_eq_mut(pos_filter)./sharp_vec_eq_orig(pos_filter);

%% Let's see what the predicted sharpness shift looks like for different
m_ind = 26;
fold_vec_eq = sharpness_vec_eq(:,m_ind) ./ sharpness_vec_eq(:,1);
fold_vec_neq = sharpness_vec_neq(:,m_ind) ./ sharpness_vec_neq(:,1);

fold_r_vec_eq = rate_vec_eq(:,m_ind) ./ rate_vec_eq(:,1);
fold_r_vec_neq = rate_vec_neq(:,m_ind) ./ rate_vec_neq(:,1);

% close all
% figure;
% hold on
% plot(cw_vec,fold_vec_eq)
% plot(cw_vec,fold_vec_neq)
% set(gca,'yscale','log')
% % xlim([2 6])