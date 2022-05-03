% script to track information rate as a function of cw

clear 
close all
addpath(genpath('../utilities/'))

% set basic parameters
nStates = 6;
paramBounds = repmat([-10 ; 6],1,3*nStates-4); % constrain transition rate magnitude
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
cw_vec = 0:0.5:8;

% cw = 1e5; %
master_struct = struct;
for c = 1:length(cw_vec)
    tic
    [sim_info_neq, sim_struct_neq] = param_sweep_multi_v3([ir_index sharpness_index],functionPath,sweep_options{:},...
                                              'half_max_flag',false,'cw',10.^cw_vec(c),...
                                              'equilibrium_flag',false,'specFactor',alpha_factor);

    [sim_info_eq, sim_struct_eq] = param_sweep_multi_v3([ir_index sharpness_index],functionPath,sweep_options{:},...
                                              'half_max_flag',false,'cw',10.^cw_vec(c),...
                                              'equilibrium_flag',true,'specFactor',alpha_factor);    
                                            
    toc
                                            
    % extract key pieces of info
    % eq
    metric_array_eq = vertcat(sim_struct_eq.metric_array);
    rate_array_eq = vertcat(sim_struct_eq.rate_array);
    
%     cw_vec_eq = metric_array_eq(:,cw_index);
    sharp_vec_eq = metric_array_eq(:,sharpness_index);
    ir_vec_eq = metric_array_eq(:,ir_index);
    ir_99_eq = prctile(ir_vec_eq,99);
    [max_sharp_eq,max_i_eq] = nanmax(sharp_vec_eq.*(1*ir_vec_eq>=ir_99_eq));
    master_struct(c).metric_array_eq = metric_array_eq;
    master_struct(c).rate_array_eq = rate_array_eq;
    master_struct(c).opt_rates_eq = rate_array_eq(max_i_eq,:);
    master_struct(c).opt_ir_eq = ir_vec_eq(max_i_eq);
    master_struct(c).opt_s_eq = max_sharp_eq;
    master_struct(c).opt_i_eq = max_i_eq;
    
    % neq (global)
    metric_array_neq = vertcat(sim_struct_neq.metric_array);
    sharp_vec_neq = metric_array_neq(:,sharpness_index);
    ir_vec_neq = metric_array_neq(:,ir_index);
    rate_array_neq = vertcat(sim_struct_neq.rate_array);
    
    ir_99_neq = prctile(ir_vec_neq,99);
    [max_sharp_neq,max_i_neq] = nanmax(sharp_vec_neq.*(1*ir_vec_neq>=ir_99_neq));
    master_struct(c).metric_array_neq = metric_array_neq;
    master_struct(c).rate_array_neq = rate_array_neq;
    master_struct(c).opt_rates_neq = rate_array_neq(max_i_neq,:);
    master_struct(c).opt_ir_neq = ir_vec_neq(max_i_neq);
    master_struct(c).opt_s_neq = sharp_vec_neq(max_i_neq);
    master_struct(c).opt_i_neq = max_i_neq;
                                            
    
end
% [sim_info_s0, sim_struct_s0] = param_sweep_multi_v3([cw_index sharp_right_norm_index],functionPath,sweep_options{:},...
%                                           'half_max_flag',false,...'cw',cw,...
%                                           'equilibrium_flag',false,'specFactor',alpha_factor);                                                                                
% 
% [sim_info_eq, sim_struct_eq] = param_sweep_multi_v3([cw_index ir_index],functionPath,sweep_options{:},...
%                                           'half_max_flag',false,...'cw',cw,...
%                                           'equilibrium_flag',true,'specFactor',alpha_factor);                                        
% toc     
%%
close all
figure;
plot([master_struct.opt_ir_neq]);
hold on
plot([master_struct.opt_ir_eq]);

figure;
plot([master_struct.opt_s_neq]);
hold on
plot([master_struct.opt_s_eq]);
%% Simulate perturbations to binding site of differing levels of severity
% add functions to path
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));

% % set basic parameters
% n_perturbations = 1e2;
% m_factor = logspace(0,log10(alpha_factor),n_perturbations);
% alpha_adjusted = alpha_factor ./ m_factor;
% c_vec = logspace(-3,6,1e4);
% c_index = find(c_vec==1);
% 
% 
% % Identify best performer in and out of equilibrium and
% 
% % eq
% metric_array_eq = vertcat(sim_struct_eq.metric_array);
% cw_vec_eq = metric_array_eq(:,cw_index);
% sharp_vec_eq = metric_array_eq(:,sharpness_index);
% ir_vec_eq = metric_array_eq(:,ir_index);
% rate_array_eq = vertcat(sim_struct_eq.rate_array);
% 
% % neq (global)
% metric_array_neq = vertcat(sim_struct_neq.metric_array);
% cw_vec_neq = metric_array_neq(:,cw_index);
% sharp_vec_neq = metric_array_neq(:,sharpness_index);
% ir_vec_neq = metric_array_neq(:,ir_index);
% rate_array_neq = vertcat(sim_struct_neq.rate_array);
% 
% % s0
% metric_array_s0 = vertcat(sim_struct_s0.metric_array);
% cw_vec_s0 = metric_array_s0(:,cw_index);
% sharp_vec_s0 = metric_array_s0(:,sharpness_index);
% ir_vec_s0 = metric_array_s0(:,ir_index);
% s0_vec_s0 = metric_array_s0(:,sharp_right_norm_index);
% f0_vec_s0 = metric_array_s0(:,spec_index);
% rate_array_s0 = vertcat(sim_struct_s0.rate_array);
% 
% % f0
% metric_array_f0 = vertcat(sim_struct_f0.metric_array);
% cw_vec_f0 = metric_array_f0(:,cw_index);
% sharp_vec_f0 = metric_array_f0(:,sharpness_index);
% ir_vec_f0 = metric_array_f0(:,ir_index);
% s0_vec_f0 = metric_array_f0(:,sharp_right_norm_index);
% f0_vec_f0 = metric_array_f0(:,spec_index);
% rate_array_f0 = vertcat(sim_struct_f0.rate_array);
% 
% % create filters
% s0_indices = f0_vec_s0<=1 & s0_vec_s0 > 1.9;
% f0_indices = s0_vec_f0<=1 & f0_vec_f0 > 1.9;
% 
% max_sharp_i_eq = NaN(1,length(cw_vec)-1);
% max_sharp_i_neq = NaN(1,length(cw_vec)-1);
% max_sharp_v_eq = NaN(1,length(cw_vec)-1);
% max_sharp_v_neq = NaN(1,length(cw_vec)-1);
% 
% max_sharp_i_s0 = NaN(1,length(cw_vec)-1);
% max_sharp_i_f0 = NaN(1,length(cw_vec)-1);
% max_sharp_v_s0 = NaN(1,length(cw_vec)-1);
% max_sharp_v_f0 = NaN(1,length(cw_vec)-1);
% 
% max_ir_v_eq = NaN(1,length(cw_vec)-1);
% max_ir_v_neq = NaN(1,length(cw_vec)-1);
% max_ir_v_s0 = NaN(1,length(cw_vec)-1);
% max_ir_v_f0 = NaN(1,length(cw_vec)-1);
% 
% for c = 1:length(cw_vec)-1
%     % obtain filters
%     cw_filter_eq = cw_vec_eq>=cw_vec(c) & cw_vec_eq<cw_vec(c+1) & sharp_vec_eq>=0;
%     cw_filter_neq = cw_vec_neq>=cw_vec(c) & cw_vec_neq<cw_vec(c+1) & sharp_vec_neq>=0;
%     cw_filter_s0 = cw_vec_s0>=cw_vec(c) & cw_vec_s0<cw_vec(c+1) & s0_indices & sharp_vec_s0>=0;
%     cw_filter_f0 = cw_vec_f0>=cw_vec(c) & cw_vec_f0<cw_vec(c+1) & f0_indices & sharp_vec_f0>=0;
%     
%     % find max value
%     [max_sharp_v_eq(c),max_sharp_i_eq(c)] = nanmax((1*cw_filter_eq).*ir_vec_eq);
%     [max_sharp_v_neq(c),max_sharp_i_neq(c)] = nanmax((1*cw_filter_neq).*ir_vec_neq);
%     [max_sharp_v_s0(c),max_sharp_i_s0(c)] = nanmax((1*cw_filter_s0).*ir_vec_s0);
%     [max_sharp_v_f0(c),max_sharp_i_f0(c)] = nanmax((1*cw_filter_f0).*ir_vec_f0);
%     
%     max_sharp_v_eq(c) = sharp_vec_eq(max_sharp_i_eq(c));
%     max_sharp_v_neq(c) = sharp_vec_neq(max_sharp_i_neq(c));
%     max_sharp_v_s0(c) = sharp_vec_s0(max_sharp_i_s0(c));
%     max_sharp_v_f0(c) = sharp_vec_f0(max_sharp_i_f0(c));
% end

%% Now simulate perturbations and observe shifts in sharpness and Kd


sharp_rate_array_eq = repmat(permute(vertcat(master_struct.opt_rates_eq),[3 2 1]), length(m_factor),1);
sharp_rate_array_eq(:,sim_info_eq.unbindingFlags==1,:) = sharp_rate_array_eq(:,sim_info_eq.unbindingFlags==1,:) .* m_factor';
sharp_rate_array_eq(:,sim_info_eq.b_index,:) = repmat(alpha_adjusted',1,1,length(cw_vec));

sharp_rate_array_neq = repmat(permute(vertcat(master_struct.opt_rates_neq),[3 2 1]), length(m_factor),1);
sharp_rate_array_neq(:,sim_info_neq.unbindingFlags==1,:) = sharp_rate_array_neq(:,sim_info_eq.unbindingFlags==1,:) .* m_factor';
sharp_rate_array_neq(:,sim_info_neq.b_index,:) = repmat(alpha_adjusted',1,1,length(cw_vec));


%%
% calculate half-max point and sharpness 
% kd_vec_eq = NaN(1,n_perturbations);
% kd_vec_neq = NaN(1,n_perturbations);
sharpness_vec_eq = NaN(n_perturbations,length(cw_vec)-1);
sharpness_vec_neq = NaN(n_perturbations,length(cw_vec)-1);

for p = 1:length(cw_vec)
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
    
end

%% Let's get a broader view of what is possible at equilibrium
[sim_info_eq, sim_struct_eq] = param_sweep_multi_v3([cw_index rate_index],functionPath,sweep_options{:},...
                                              'half_max_flag',false,...
                                              'equilibrium_flag',true,'specFactor',alpha_factor); 

%%                                            
rates_eq = vertcat(sim_struct_eq.rate_array);
metrics_eq = vertcat(sim_struct_eq.metric_array);
cw_vec_eq = metrics_eq(:,cw_index);

m_factor_index = 51; 
rates_eq_mutated = rates_eq;
rates_eq_mutated(:,sim_info_eq.unbindingFlags==1) = rates_eq_mutated(:,sim_info_eq.unbindingFlags==1) * m_factor(m_factor_index);
rates_eq_mutated(sim_info_eq.b_index,:) = alpha_adjusted(m_factor_index);

paramCellEqMut = mat2cell(rates_eq_mutated,size(rates_eq_mutated,1),ones(1,size(rates_eq_mutated,2)));
sharp_vec_eq_orig = metrics_eq(:,sharpness_index);
sharp_vec_eq_mut = sharpnessFunction(paramCellEqMut{:});
pos_filter = sharp_vec_eq_mut>=0 & sharp_vec_eq_orig >= 0;

fold_change_vec = sharp_vec_eq_mut(pos_filter)./sharp_vec_eq_orig(pos_filter);

%% Let's see what the predicted sharpness shift looks like for different

cw_plot_values = log10([1 alpha_factor alpha_factor^2 alpha_factor^3]);
cw_plot_indices = discretize(cw_plot_values,cw_vec);
fold_vec_eq = sharpness_vec_eq(51,:) ./ sharpness_vec_eq(1,:);
fold_vec_neq = sharpness_vec_neq(51,:) ./ sharpness_vec_neq(1,:);

close all
figure;
hold on
plot(cw_vec,fold_vec_eq)
plot(cw_vec,fold_vec_neq)
set(gca,'yscale','log')
% xlim([2 6])