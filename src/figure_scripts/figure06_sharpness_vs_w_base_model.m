% script to track observed sharpness as a function of CW
clear 
close all
addpath(genpath('../utilities/'))

% set basic parameters
nStates = 6;
paramBounds = repmat([-5 ; 5],1,11); % constrain transition rate magnitude
paramBounds(1,8) = 0; % ensures activation
paramBounds(2,9) = 0; % ensures activation
[~,~,metric_names] = calculateMetricsSym_v2([]);

% specify function path
% functionPath = ['../utilities/metricFunctions/n' num2str(nStates) '_OR/'];
functionPath = '../utilities/metricFunctions/symbolic/n006_s01_ns01_g01';

% make sure we're linked to the appropriate function subfolder% make sure we're linked to the appropriate function subfolder
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));

% define save path

DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\manuscript\';
% DropboxFolder = 'S:\Nick\Dropbox\Nonequilibrium\Nick\manuscript\';

FigPath = [DropboxFolder 'experimental_signatures' filesep];
mkdir(FigPath);         

% get index of useful metrics
spec_index = find(strcmp(metric_names,'Specificity'));
sharp_right_norm_index = find(strcmp(metric_names,'SharpnessRightNorm'));
ir_index = find(strcmp(metric_names,'DecisionRateNorm'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
precision_index = find(strcmp(metric_names,'Precision'));
rate_index = find(strcmp(metric_names,'Production Rate'));
precision_right_index = find(strcmp(metric_names,'PrecisionRight'));
cw_index = find(strcmp(metric_names,'CW'));

% set sim options
sweep_options = {'n_seeds',5,'n_iters_max',50,'n_sim',10,'nStates',nStates};

% calculate sensitivity bound
alpha_factor = 100;
% seq = 1/4;

% specify plot options 
markerAlpha = 0.5; % marker transparency
markerSize = 75; % marker size
% n_plot = 3e3;
% bit_factor = log2(exp(1));
%% 
% observed sharpness vs cw
tic
[sim_info_neq, sim_struct_neq] = param_sweep_multi_v3([cw_index sharpness_index],functionPath,sweep_options{:},...
                                          'half_max_flag',false,...
                                          'equilibrium_flag',false,'specFactor',alpha_factor,'paramBounds',paramBounds);

[sim_info_eq, sim_struct_eq] = param_sweep_multi_v3([cw_index sharpness_index],functionPath,sweep_options{:},...
                                          'half_max_flag',false,...
                                          'equilibrium_flag',true,'specFactor',alpha_factor,'paramBounds',paramBounds);                                        
toc     


% seq = 1/4;
cr = 1;

%% Generate bound predictions
rd = [190 30 45]/255;

cw_axis = logspace(0,5);
s_bound_fun = @(f) (f./(f+cw_axis)).*(1 + (alpha_factor^2-f)./(f.*alpha_factor + cw_axis.*(f-alpha_factor)));
s0_max_pd = s_bound_fun(alpha_factor);
f0_max_pd = s_bound_fun(alpha_factor^2);
eq_max_pd = (alpha_factor./(alpha_factor+cw_axis));

% extract vectors to plot (neq)
n_plot = 1e4;
results_array_neq = vertcat(sim_struct_neq.metric_array);
cw_scatter_vec_neq = results_array_neq(:,cw_index);
s_scatter_vec_neq = results_array_neq(:,sharpness_index);
plot_options_cw_neq = find(s_scatter_vec_neq>=0);
plot_indices_cw_neq = randsample(plot_options_cw_neq,min([n_plot,length(plot_options_cw_neq)]),false);

% extract vectors to plot (eq)
results_array_eq = vertcat(sim_struct_eq.metric_array);
cw_scatter_vec_eq = results_array_eq(:,cw_index);
s_scatter_vec_eq = results_array_eq(:,sharpness_index);
plot_options_cw_eq = find(s_scatter_vec_eq>=0);
plot_indices_cw_eq = randsample(plot_options_cw_eq,min([n_plot,length(plot_options_cw_eq)]),false);

%%%%%%%%%%%%%%%%%
% generate predicted curves for the two motifs
cw_vec = logspace(nanmin(cw_scatter_vec_neq),nanmax(cw_scatter_vec_neq),101);
% cw = 1./cw_vec;
s_bound_fun = @(f0) ((-1)+alpha_factor).^(-1).*(alpha_factor.^2+((-2)+alpha_factor).*f0).*(cw_vec+cr.*f0).^(-1);
                
% s_bound_fun = @(f0) (f0./cw_vec) ./ (1 + f0./cw_vec) .* seq .* (alpha_factor-1) ./ alpha_factor .* (alpha_factor + f0) ./ (f0- 1); 
% s_bound_s0 = seq*s_bound_fun(alpha_factor);
% s_bound_f0 = seq*s_bound_fun(alpha_factor^2);
% s_bound_eq = seq*alpha_factor ./ (alpha_factor + cw_vec);
% cw_axis = (cw_vec) ;

%%%%%%%%%%%%%%%%%%%%%%%%
% make figure
close all
motif_fig = figure;
hold on
cmap = brewermap([],'Set2');

% plot numerical results
sneq = scatter(10.^cw_scatter_vec_neq(plot_indices_cw_neq),s_scatter_vec_neq(plot_indices_cw_neq),...
      markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha*0.5, 'MarkerFaceColor',cmap(2,:));
    
seq = scatter(10.^cw_scatter_vec_eq(plot_indices_cw_eq),s_scatter_vec_eq(plot_indices_cw_eq),...
      markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha*0.5, 'MarkerFaceColor',cmap(3,:));    

% plot(cw_axis,s0_max_pd,'-.','Color',cmap(5,:),'LineWidth',3)
% plot(cw_axis,f0_max_pd,'-.','Color',cmap(4,:),'LineWidth',3)


% s0-optimized prediction
% p1 = scatter(cw_axis,s_bound_f0,'^','MarkerFaceColor',cmap(5,:),'MarkerEdgeColor','k');
% f0-optimized prediction
% p2 = scatter(cw_axis,s_bound_s0,'s','MarkerFaceColor',cmap(4,:),'MarkerEdgeColor','k');
% % equilibrium bound
% p3 = scatter(cw_axis,s_bound_eq,'o','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k');


set(gca,'xscale','log')
xlabel('relative non-cognate factor concentration (w / c)');
ylabel('normalized sharpness (S)')

set(gca,'xtick',[1 10^1 10^2 10^3 10^4 10^5]);%,'xticklabels',{'\alpha^0','\alpha^1','\alpha^2','\alpha^3','\alpha^4'})
% legend([p2 p1 s1],'sharpness optimized','fidelity optimized','global optimum')

xlim([1e0 10^5])
ylim([0 2])
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

motif_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

% set(p3,'Visible','off')
% saveas(motif_fig,[FigPath 'effective_sharpness_vs_cw_scatter_s0_f0.png'])
saveas(motif_fig,[FigPath 'sharpness_vs_cw_scatter.png'])
% saveas(motif_fig,[FigPath 'sharpness_vs_cw_scatter.pdf'])

delete(seq)
delete(sneq)


plot(cw_axis,max([f0_max_pd ; s0_max_pd]),'-.','Color',rd,'LineWidth',3)
plot(cw_axis,eq_max_pd,'-.','Color',brighten(cmap(3,:),-0.6),'LineWidth',3)
grid on
% set(p3,'Visible','on')
saveas(motif_fig,[FigPath 'sharpness_vs_cw_scatter_lines.pdf'])


%% Plot naive and actual sharpness bounds at reasonable cw level
cw_true = 1e3;
cw_naive = 1e-2; % doesn't matter so long as leq 1
sweep_options = {'n_seeds',5,'n_iters_max',50,'n_sim',50,'nStates',nStates};
% "true" simulations
tic
[sim_info_neq, sim_struct_neq] = param_sweep_multi_v3([sharpness_index precision_index],functionPath,sweep_options{:},...
                                          'half_max_flag',false,'cw',cw_true,...
                                          'equilibrium_flag',false,'specFactor',alpha_factor,'paramBounds',paramBounds);    

[sim_info_eq, sim_struct_eq] = param_sweep_multi_v3([sharpness_index precision_index],functionPath,sweep_options{:},...
                                          'half_max_flag',false,'cw',cw_true,...
                                          'equilibrium_flag',true,'specFactor',alpha_factor,'paramBounds',paramBounds);                                    
toc     


% now ask what happens when we neglect to account for CW
[sim_info_eq_naive, sim_struct_eq_naive] = param_sweep_multi_v3([sharpness_index precision_index],functionPath,sweep_options{:},...
                                          'half_max_flag',false,'cw',cw_naive,...
                                          'equilibrium_flag',true,'specFactor',alpha_factor,'paramBounds',paramBounds);        

close all
%%
rate_bounds = [0.49 0.51];%-0.25;

% Calculate titration curve for each scenario
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));

% extract metric arrays
metric_array_neq = sim_struct_neq.metric_array;
rate_vec_neq = metric_array_neq(:,rate_index);
hm_filter_neq = rate_vec_neq>=rate_bounds(1)&rate_vec_neq<=rate_bounds(2);
metric_array_eq = sim_struct_eq.metric_array;
rate_vec_eq = metric_array_eq(:,rate_index);
hm_filter_eq = rate_vec_eq>=rate_bounds(1)&rate_vec_eq<=rate_bounds(2);
metric_array_eq_naive = sim_struct_eq_naive.metric_array;
rate_vec_eq_naive = metric_array_eq_naive(:,rate_index);
hm_filter_eq_naive = rate_vec_eq_naive>=rate_bounds(1)&rate_vec_eq_naive<=rate_bounds(2);

%% rate arrays
rate_array_neq = vertcat(sim_struct_neq.rate_array);
rate_array_eq = vertcat(sim_struct_eq.rate_array);
rate_array_eq_naive = vertcat(sim_struct_eq_naive.rate_array);

% identify sharpest networks
[max_s_neq,max_s_i_neq] = nanmax(metric_array_neq(:,sharpness_index).*hm_filter_neq);
[max_s_eq,max_s_i_eq] = nanmax(metric_array_eq(:,sharpness_index).*hm_filter_eq);
[max_s_eq_naive,max_s_i_eq_naive] = nanmax(metric_array_eq_naive(:,sharpness_index).*hm_filter_eq_naive);

% c vector
cr_vec = logspace(-2,2,1e3);

n_plot = 50;
ir_vec_neq = metric_array_neq(:,ir_index);
s_vec_neq = metric_array_neq(:,sharpness_index);
ir_98 = prctile(ir_vec_neq(hm_filter_neq),99);
rng(123);
plot_indices_neq = randsample(find(ir_vec_neq>ir_98&s_vec_neq>=0&hm_filter_neq),n_plot,false);

opt_rates_neq = NaN(length(cr_vec),n_plot);

for n = 1:n_plot
    temp_array = repmat(rate_array_neq(plot_indices_neq(n),:),length(cr_vec),1);
    temp_array(:,sim_info_neq.cr_index) = cr_vec;
    paramCellNeq = mat2cell(temp_array,size(temp_array,1),ones(1,size(temp_array,2)));
    opt_rates_neq(:,n) = productionRateFunction(paramCellNeq{:});
end


% calculate production rates
eq_rate_array = NaN(length(cr_vec),2);

best_eq_rate_array = repmat(rate_array_eq(max_s_i_eq,:),length(cr_vec),1);
best_eq_rate_array(:,sim_info_eq.cr_index) = cr_vec;
paramCellEq = mat2cell(best_eq_rate_array,size(best_eq_rate_array,1),ones(1,size(best_eq_rate_array,2)));
eq_rate_array(:,1) = productionRateFunction(paramCellEq{:});

best_eq_naive_rate_array = repmat(rate_array_eq_naive(max_s_i_eq_naive,:),length(cr_vec),1);
best_eq_naive_rate_array(:,sim_info_eq_naive.cr_index) = cr_vec;
paramCellEqNaive = mat2cell(best_eq_naive_rate_array,size(best_eq_naive_rate_array,1),ones(1,size(best_eq_naive_rate_array,2)));
eq_rate_array(:,2) = productionRateFunction(paramCellEqNaive{:});

% numbers for paper
H_neq_low = prctile(metric_array_neq(plot_indices_neq,sharpness_index),25)
H_neq_high = prctile(metric_array_neq(plot_indices_neq,sharpness_index),75)
H_eq_max = max(metric_array_eq(:,sharpness_index))
%% Make Figures

titration_fig = figure;
hold on
cmap = brewermap([],'Set2');

sneq = plot(cr_vec, opt_rates_neq,'Color',[cmap(2,:) .5],'LineWidth',1);
seq = plot(cr_vec, eq_rate_array(:,1),'Color',cmap(3,:),'LineWidth',3);
seq_naive = plot(cr_vec, eq_rate_array(:,2),'--','Color',cmap(3,:),'LineWidth',3);
xlim([0.1 10])    
yl = get(gca,'ylim');
% set(gca,'yscale','log')
% yl(end) = 10;
plot(repmat([0.95 1.05],2,1),repmat(yl',1,2),'-.','Color','k','LineWidth',1.5)

% xlim([alpha_factor alpha_factor^2])
set(gca,'xscale','log')    
xlabel('activator concentration (C)');
ylabel('production rate (R)')
%legend([sneq(1) seq_naive seq],'non-equilibrium','equilibrium (naive)','equilibrium (actual)','Location','southeast')

grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

titration_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(titration_fig,[FigPath 'sharpness_naive_vs_reality.png'])
saveas(titration_fig,[FigPath 'sharpness_naive_vs_reality.pdf'])

% %% %%%%%%%%%%%%%%%%%%%%%%%%%%% info vs cw %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % randomly draw ~100 top performing non-eq networks
% close all
% 
% % calculate actual and niave equilibrium boundaries
% eq_array = [metric_array_eq(:,sharpness_index), exp(metric_array_eq(:,precision_index))];
% eq_array = eq_array(all(~isnan(eq_array),2)&eq_array(:,1)>=0,:);
% 
% naive_eq_array = [metric_array_eq_naive(:,sharpness_index), exp(metric_array_eq_naive(:,precision_index))];
% naive_eq_array = naive_eq_array(all(~isnan(naive_eq_array),2)&naive_eq_array(:,1)>=0,:);
% 
% % eq_bound_points = boundary(eq_array,1);
% % naive_eq_bound_points = boundary(naive_eq_array,1);
% 
% eq_max = nanmax(eq_array);
% mval = 1e-6;
% eq_bound_points = [mval,eq_max(2) ; eq_max ; eq_max(1) mval;mval mval];
% naive_eq_max = nanmax(naive_eq_array);
% naive_eq_bound_points = [mval,naive_eq_max(2) ; naive_eq_max; naive_eq_max(1) mval; mval mval];
% 
% %% Make plot
% phase_space_fig = figure;
% hold on
% fill([naive_eq_bound_points(:,1)' fliplr(naive_eq_bound_points(:,1)')],...
%                 [naive_eq_bound_points(:,2)' naive_eq_bound_points(:,2)'],cmap(3,:),'LineWidth',2,'LineStyle','--','FaceAlpha',0.25)
% fill([eq_bound_points(:,1)' fliplr(eq_bound_points(:,1)')],...
%                 [eq_bound_points(:,2)' eq_bound_points(:,2)'],cmap(3,:),'LineWidth',2,'FaceAlpha',0.75)
% 
% scatter(metric_array_neq(plot_indices_neq,sharpness_index),metric_array_neq(plot_indices_neq,precision_index),...
%                         'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k');
% grid on
% set(gca,'FontSize',14)
% set(gca,'Color',[228,221,209]/255) 
% ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.XAxis(1).Color = 'k';
% xlim([1e-4 1])
% ylim([1e-4 1e2])
% phase_space_fig.InvertHardcopy = 'off';
% set(gcf,'color','w');
% 
% xlabel('observed sharpness (s)')
% ylabel('precision (1/\sigma^2)')
% 
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% saveas(phase_space_fig,[FigPath 'phase_space_naive_vs_real.png'])
% saveas(phase_space_fig,[FigPath 'phase_space_naive_vs_real.pdf'])