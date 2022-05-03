% script to generate data and figures for fidelity v sharpness parameter
% space exploration

clear 
close all
addpath('../utilities/')
FigPath = '../../fig/bivariate_param_sweeps/';
mkdir(FigPath);
[~, metric_names] = calculateMetricsV2([],[],[]);

% call edge sampler for fidelity vs. absolute range
metric_indices = [5,8];
[sim_struct_eq, fpath_eq] = param_sweep(metric_indices,1);
[sim_struct_neq, fpath_neq] = param_sweep(metric_indices,0);

%% make plots
type_names = {'burst frequency', 'burst duration'};
% plot metric dispersion
mSize = 25;
close all
cmap = brewermap(9,'Set2');
color_indices = [2,3,5,6 ; 1,4,9,10];
n_plot = 2000;
for i = 1%:numel(sim_struct_eq)        
    close all
    metric_array_eq = sim_struct_eq(i).metric_array;
    metric_array_neq = sim_struct_neq(i).metric_array;   
    
    
    metric_fig = figure;
    hold on                   
    % extract vectors
    mvec1_eq = reshape(metric_array_eq(:,metric_indices(1),:,:),[],1);
    mvec2_eq = 1./reshape(metric_array_eq(:,metric_indices(2),:,:),[],1);
    mvec1_neq = reshape(metric_array_neq(:,metric_indices(1),:,:),[],1);
    mvec2_neq = 1./reshape(metric_array_neq(:,metric_indices(2),:,:),[],1);
    index_vec = 1:numel(mvec1_eq);    
    plot_indices = randsample(index_vec,n_plot,false);
    % plot eq and non-eq scatters       
    b_neq = convhull(mvec1_neq,mvec2_neq);
    b_eq = convhull(mvec1_eq,mvec2_eq);
    scatter(mvec1_neq(plot_indices),mvec2_neq(plot_indices),mSize,'MarkerFaceColor',cmap(color_indices(i,1),:),'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2);        
    scatter(mvec1_eq(plot_indices),mvec2_eq(plot_indices),mSize,'MarkerFaceColor',cmap(color_indices(i,2),:),'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2);        
    pneq = plot(mvec1_neq(b_neq),mvec2_neq(b_neq),'Color',cmap(color_indices(i,1),:),'LineWidth',3);    
    peq = plot(mvec1_eq(b_eq),mvec2_eq(b_eq),'Color',cmap(color_indices(i,2),:),'LineWidth',3);

    grid on
    box on    
    xlabel('normalized fold-change')
    ylabel('normalized mRNA production rate')
    legend([pneq peq], 'non-equilibrium','equilibrium', 'Location','northeast','Fontsize',12) 
    saveas(metric_fig,[FigPath type_names{i} '_' metric_names{metric_indices(1)} '_' metric_names{metric_indices(2)} '.png'])
    saveas(metric_fig,[FigPath type_names{i} '_' metric_names{metric_indices(1)} '_' metric_names{metric_indices(2)} '.pdf'])
end

% close all