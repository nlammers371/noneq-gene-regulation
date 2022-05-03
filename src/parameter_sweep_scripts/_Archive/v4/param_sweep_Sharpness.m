% script to generate data and figures for fidelity v sharpness parameter
% space exploration
clear 
close all
addpath('../utilities/')
FigPath = '../../fig/bivariate_param_sweeps/';
mkdir(FigPath);
[~, metric_names] = calculateMetricsV2([],[],[]);
% call edge sampler for production rate vs. precision
metric_indices_info = [7,8];
[sim_struct_info_eq, ~] = param_sweep(metric_indices_info,1,'n_seeds',10);
[sim_struct_info_neq, ~] = param_sweep(metric_indices_info,0,'n_seeds',10);

% call edge sampler for sharpness vs. production rate
metric_indices_pd = [4,7];
[sim_struct_pd_eq, ~] = param_sweep(metric_indices_pd,1);
[sim_struct_pd_neq, ~] = param_sweep(metric_indices_pd,0);

% call edge sampler for sharpness vs. production rate
metric_indices_noise = [4,6];
[sim_struct_ns_eq, ~] = param_sweep(metric_indices_noise,1,'n_seeds',10);
[sim_struct_ns_neq, ~] = param_sweep(metric_indices_noise,0,'n_seeds',10);

%% make plots
type_names = {'burst frequency', 'burst duration'};
% plot metric dispersion
mSize = 25;
close all
cmap = brewermap(12,'Paired');
color_indices = [1,2,5,6 ; 3,4,9,10];
n_plot = 5000;
for i = 1:numel(sim_struct_pd_eq)        
    close all    
    
    %%%%%%%%%% production rate vs. sharpness
    metric_array_pd_eq = sim_struct_pd_eq(i).metric_array;
    metric_array_pd_neq = sim_struct_pd_neq(i).metric_array;        
    
    mvec1_pd_eq = reshape(metric_array_pd_eq(:,metric_indices_pd(1),:,:),[],1);
    mvec2_pd_eq = reshape(metric_array_pd_eq(:,metric_indices_pd(2),:,:),[],1);
    mvec1_pd_neq = reshape(metric_array_pd_neq(:,metric_indices_pd(1),:,:),[],1);
    mvec2_pd_neq = reshape(metric_array_pd_neq(:,metric_indices_pd(2),:,:),[],1);
    b_neq = convhull(mvec1_pd_neq,mvec2_pd_neq );
    b_eq = convhull(mvec1_pd_eq,mvec2_pd_eq);
    % draw indices to plot
    index_vec = 1:numel(mvec1_pd_eq);    
    n_plot = 5000;
    plot_indices = randsample(index_vec,n_plot,false);
    % plot eq and non-eq scatters       
    
    pd_sharpness = figure;
    hold on    
    scatter(mvec1_pd_neq(plot_indices),mvec2_pd_neq(plot_indices),mSize,'MarkerFaceColor',cmap(color_indices(i,1),:),'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2);            
    scatter(mvec1_pd_eq(plot_indices),mvec2_pd_eq(plot_indices),mSize,'MarkerFaceColor',cmap(color_indices(i,3),:),'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2);        
    pneq = plot(mvec1_pd_neq(b_neq),mvec2_pd_neq(b_neq),'Color',cmap(color_indices(i,2),:),'LineWidth',3);
    peq = plot(mvec1_pd_eq(b_eq),mvec2_pd_eq(b_eq),'Color',cmap(color_indices(i,4),:),'LineWidth',3);
    grid on
    box on    
    xlabel('sharpness (units of 1/c)')
    ylabel('normalized mRNA production rate')
    set(gca,'FontSize',14)
    legend([pneq peq], 'non-equilibrium','equilibrium', 'Location','northeast','Fontsize',12) 
    saveas(pd_sharpness,[FigPath type_names{i} '_' metric_names{metric_indices_pd(1)} '_' metric_names{metric_indices_pd(2)} '.png'])
    saveas(pd_sharpness,[FigPath type_names{i} '_' metric_names{metric_indices_pd(1)} '_' metric_names{metric_indices_pd(2)} '.pdf'])    
    
    %%%%%%%%%% production rate vs. precision
    metric_array_info_eq = sim_struct_info_eq(i).metric_array;
    metric_array_info_neq = sim_struct_info_neq(i).metric_array;        
    
    mvec1_info_eq = reshape(metric_array_info_eq(:,metric_indices_info(1),:,:),[],1);
    mvec2_info_eq = reshape(metric_array_info_eq(:,metric_indices_info(2),:,:),[],1);
    mvec1_info_neq = reshape(metric_array_info_neq(:,metric_indices_info(1),:,:),[],1);
    mvec2_info_neq = reshape(metric_array_info_neq(:,metric_indices_info(2),:,:),[],1);
    b_neq = convhull(mvec1_info_neq,mvec2_info_neq);
    b_eq = convhull(mvec1_info_eq,mvec2_info_eq);
    index_vec = 1:numel(mvec1_info_eq);    
    n_plot = 5000;
    plot_indices = randsample(index_vec,n_plot,false);
    pd_info = figure;
    hold on    
    scatter(mvec1_info_neq(plot_indices),mvec2_info_neq(plot_indices),mSize,'MarkerFaceColor',cmap(color_indices(i,1),:),'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2);            
    scatter(mvec1_info_eq(plot_indices),mvec2_info_eq(plot_indices),mSize,'MarkerFaceColor',cmap(color_indices(i,3),:),'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2);        
    pneq = plot(mvec1_info_neq(b_neq),mvec2_info_neq(b_neq),'Color',cmap(color_indices(i,2),:),'LineWidth',3);
    peq = plot(mvec1_info_eq(b_eq),mvec2_info_eq(b_eq),'Color',cmap(color_indices(i,4),:),'LineWidth',3);
    grid on
    box on    
    ylabel('precision (1/\sigma_c)')
    xlabel('normalized mRNA production rate')
    set(gca,'FontSize',14)
    legend([pneq peq], 'non-equilibrium','equilibrium', 'Location','northeast','Fontsize',12) 
    saveas(pd_info,[FigPath type_names{i} '_' metric_names{metric_indices_info(1)} '_' metric_names{metric_indices_info(2)} '.png'])
    saveas(pd_info,[FigPath type_names{i} '_' metric_names{metric_indices_info(1)} '_' metric_names{metric_indices_info(2)} '.pdf'])    
    
    %%%%%%%%%% production rate vs. precision
    metric_array_ns_eq = sim_struct_ns_eq(i).metric_array;
    metric_array_ns_neq = sim_struct_ns_neq(i).metric_array;        
    
    mvec1_ns_eq = reshape(metric_array_ns_eq(:,metric_indices_noise(1),:,:),[],1);
    mvec2_ns_eq = reshape(metric_array_ns_eq(:,metric_indices_noise(2),:,:),[],1);
    mvec1_ns_neq = reshape(metric_array_ns_neq(:,metric_indices_noise(1),:,:),[],1);
    mvec2_ns_neq = reshape(metric_array_ns_neq(:,metric_indices_noise(2),:,:),[],1);
    
    mvec2_ns_neq = mvec2_ns_neq(mvec1_ns_neq>.05);
    mvec1_ns_neq = mvec1_ns_neq(mvec1_ns_neq>.05);
    mvec2_ns_eq = mvec2_ns_eq(mvec1_ns_eq>.05);
    mvec1_ns_eq = mvec1_ns_eq(mvec1_ns_eq>.05);
    if i == 2
        b_neq = boundary(mvec1_ns_neq,mvec2_ns_neq,.9);
        b_eq = boundary(mvec1_ns_eq,mvec2_ns_eq,.9);
    else
        b_neq = boundary(mvec1_ns_neq,mvec2_ns_neq,.7);
        b_eq = boundary(mvec1_ns_eq,mvec2_ns_eq,.7);
    end
    index_vec = 1:min([numel(mvec1_ns_eq),numel(mvec1_ns_neq)]);    
    n_plot = 5000;
    plot_indices = randsample(index_vec,n_plot,false);
    
    pd_noise = figure;
    hold on    
    scatter(mvec1_ns_neq(plot_indices),mvec2_ns_neq(plot_indices),mSize,'MarkerFaceColor',cmap(color_indices(i,1),:),'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2);            
    scatter(mvec1_ns_eq(plot_indices),mvec2_ns_eq(plot_indices),mSize,'MarkerFaceColor',cmap(color_indices(i,3),:),'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2);        
    pneq = plot(mvec1_ns_neq(b_neq),mvec2_ns_neq(b_neq),'Color',cmap(color_indices(i,2),:),'LineWidth',3);
    peq = plot(mvec1_ns_eq(b_eq),mvec2_ns_eq(b_eq),'--','Color',cmap(color_indices(i,4),:),'LineWidth',3);
    grid on
    xlim([.05,.5])
    box on    
    ylabel('precision (1/\sigma_c)')
    xlabel('normalized mRNA production rate')
    set(gca,'FontSize',14)
    legend([pneq peq], 'non-equilibrium','equilibrium', 'Location','northeast','Fontsize',12) 
    saveas(pd_noise,[FigPath type_names{i} '_' metric_names{metric_indices_noise(1)} '_' metric_names{metric_indices_noise(2)} '.png'])
    saveas(pd_noise,[FigPath type_names{i} '_' metric_names{metric_indices_noise(1)} '_' metric_names{metric_indices_noise(2)} '.pdf'])    
end

    
