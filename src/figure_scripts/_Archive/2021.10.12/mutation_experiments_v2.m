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

FigPath = [DropboxFolder 'observed_sharpness' filesep];
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

% Load data
load([FigPath 'mutation_experiment_data'],'master_struct')

%% Simulate perturbations to binding site of differing levels of severity
% add functions to path
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));

[sim_info_neq, ~] = param_sweep_multi_v3([ir_index sharpness_index],functionPath,sweep_options{:},...
                                              'half_max_flag',false,'cw',1,...
                                              'equilibrium_flag',false,'specFactor',alpha_factor,'paramBounds',paramBounds,'n_sim',1,'n_iters_max',1);
[sim_info_eq, ~] = param_sweep_multi_v3([ir_index sharpness_index],functionPath,sweep_options{:},...
                                              'half_max_flag',false,'cw',1,...
                                              'equilibrium_flag',false,'specFactor',alpha_factor,'paramBounds',paramBounds,'n_sim',1,'n_iters_max',1);                                            
% set basic parameters
n_perturbations = 1e2;
m_factor = logspace(0,log10(alpha_factor),n_perturbations);
alpha_adjusted = alpha_factor ./ m_factor;
c_vec = logspace(-3,6,1e4);
c_index = find(c_vec==1);


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


%% Let's see what the predicted sharpness shift looks like for different

m_ind = 26;
sharpness_array_norm_eq = sharpness_vec_eq ./ sharpness_vec_eq(:,1);
sharpness_array_norm_neq = sharpness_vec_neq ./ sharpness_vec_neq(:,1);

fold_vec_eq = sharpness_vec_eq(:,m_ind) ./ sharpness_vec_eq(:,1);
fold_vec_neq = sharpness_vec_neq(:,m_ind) ./ sharpness_vec_neq(:,1);

fold_r_vec_eq = rate_vec_eq(:,m_ind) ./ rate_vec_eq(:,1);
fold_r_vec_neq = rate_vec_neq(:,m_ind) ./ rate_vec_neq(:,1) ;

cw_index = unique(cw_vec_neq);
[~,cw_plot_ind] = nanmin(abs(cw_index - 4));

fold_s_neq_plot = NaN(size(cw_index));
fold_s_eq_plot = NaN(size(cw_index));
fold_r_eq_plot = NaN(size(cw_index));
fold_r_neq_plot = NaN(size(cw_index));
for c = 1:length(cw_index)
    % fold change vs cw
    fold_s_eq_plot(c) = nanmax(fold_vec_eq(cw_vec_eq==cw_index(c)));
    fold_s_neq_plot(c) = nanmax(fold_vec_neq(cw_vec_neq==cw_index(c)));
    
    fold_r_eq_plot(c) = nanmax(fold_r_vec_eq(cw_vec_eq==cw_index(c)));
    fold_r_neq_plot(c) = nanmax(fold_r_vec_neq(cw_vec_neq==cw_index(c)));        
end    

plot_c = 4;
fold_s_cw_neq_plot = NaN(size(cw_index));
fold_s_cw_eq_plot = NaN(size(cw_index));
for m = 1:length(m_factor)
    % fold change vs alpha
    fold_s_cw_eq_plot(m) = nanmax(sharpness_array_norm_eq(cw_vec_eq==plot_c,m));
    fold_s_cw_neq_plot(m) = nanmax(sharpness_array_norm_neq(cw_vec_neq==plot_c,m));
end  

close all
% sharpness fold shift vs perturbation strength
alpha_shift_fig = figure;
cmap = brewermap([],'Set2');
hold on
% scatter(m_factor / alpha_factor,fold_s_cw_neq_plot,markerSize,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k')
% scatter(m_factor / alpha_factor,fold_s_cw_eq_plot,markerSize,'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k')
plot(m_factor / alpha_factor,fold_s_cw_neq_plot,'Color',cmap(2,:),'LineWidth',3)
plot(m_factor / alpha_factor,fold_s_cw_eq_plot,'Color',cmap(3,:),'LineWidth',3)
set(gca,'yscale','log')
set(gca,'xscale','log')
xlabel('perturbation strength (\alpha^* / \alpha)');
ylabel('fold sharpness decrease (s^*/s)')
% xlim([1e0 alpha_factor^4])
ylim([10^-4 10^0])

% set(gca,'xtick',[1 alpha_factor^1 alpha_factor^2 alpha_factor^3 alpha_factor^4]);%,'xticklabels',{'\alpha^0','\alpha^1','\alpha^2','\alpha^3','\alpha^4'})

set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
alpha_shift_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
saveas(alpha_shift_fig,[FigPath 'sharp_fold_vs_alpha.png'])
saveas(alpha_shift_fig,[FigPath 'sharp_fold_vs_alpha.pdf'])
%%
% sharpness fold shift vs cw
sharp_shift_fig = figure;
cmap = brewermap([],'Set2');
hold on
scatter(10.^cw_index,fold_s_neq_plot,markerSize,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k')
scatter(10.^cw_index,fold_s_eq_plot,markerSize,'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k')
set(gca,'yscale','log')
set(gca,'xscale','log')
xlabel('relative non-cognate factor concentration (c_w / c_r)');
ylabel('fold sharpness decrease (s^*/s)')
xlim([1e0 alpha_factor^4])
ylim([10^-1.1 10^0])

set(gca,'xtick',[1 alpha_factor^1 alpha_factor^2 alpha_factor^3 alpha_factor^4]);%,'xticklabels',{'\alpha^0','\alpha^1','\alpha^2','\alpha^3','\alpha^4'})

set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
sharp_shift_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
saveas(sharp_shift_fig,[FigPath 'sharp_fold_vs_cw.png'])
saveas(sharp_shift_fig,[FigPath 'sharp_fold_vs_cw.pdf'])

% rate shift
rate_shift_fig = figure;
cmap = brewermap([],'Set2');
hold on
scatter(10.^cw_index,fold_r_neq_plot,markerSize,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k')
scatter(10.^cw_index,fold_r_eq_plot,markerSize,'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k')
set(gca,'yscale','log')
set(gca,'xscale','log')
xlabel('relative non-cognate factor concentration (c_w / c_r)');
ylabel('fold sharpness decrease (r^*/r)')

set(gca,'xtick',[1 alpha_factor^1 alpha_factor^2 alpha_factor^3 alpha_factor^4]);%,'xticklabels',{'\alpha^0','\alpha^1','\alpha^2','\alpha^3','\alpha^4'})
% legend([p2 p1 s1],'sharpness optimized','fidelity optimized','global optimum')
xlim([1e0 alpha_factor^4])
set(gca,'FontSize',14)
ylim([1e-1 10^0.5])
set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
set(gca,'yscale','log')
rate_shift_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(rate_shift_fig,[FigPath 'rate_fold_vs_cw.png'])
saveas(rate_shift_fig,[FigPath 'rate_fold_vs_cw.pdf'])

% figure;
% hold on
% scatter(cw_vec_neq,fold_r_vec_neq)
% hold on
% scatter(cw_vec_eq,fold_r_vec_eq)
% set(gca,'yscale','log')