% script to generate data and figures for fidelity v sharpness parameter
% space exploration

clear 
close all
addpath('../utilities/')
FigPath = '../../fig/bivariate_param_sweeps/';
mkdir(FigPath);
[~, metric_names] = calculateMetricsV2([],[],[]);

% call edge sampler for fidelity vs. absolute range
metric_indices = [5,7];
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
for i = 1:numel(sim_struct_eq)        
    close all
    metric_array_eq = sim_struct_eq(i).metric_array;
    metric_array_neq = sim_struct_neq(i).metric_array;        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % first make basic proofreading fig    
    % extract vectors
    fold_in_eq = reshape(metric_array_eq(:,1,:,:),[],1);
    fold_out_eq = reshape(metric_array_eq(:,2,:,:),[],1);
    fold_in_neq = reshape(metric_array_neq(:,1,:,:),[],1);
    fold_out_neq = reshape(metric_array_neq(:,2,:,:),[],1);
    b_neq = convhull(fold_in_neq,fold_out_neq);
    b_eq = convhull(fold_in_eq,fold_out_eq);
    % draw indices to plot
    index_vec = 1:numel(fold_out_neq);    
    plot_indices = randsample(index_vec,n_plot,false);
    % calculate plot boundaries
    y_max = max(fold_out_neq(b_neq));
    x_max = max(fold_in_neq(b_neq));
    % plot eq and non-eq scatters       
    
    in_out_fold_eq = figure;
    hold on    
    scatter(fold_in_eq(plot_indices),fold_out_eq(plot_indices),mSize,'MarkerFaceColor',cmap(color_indices(i,2),:),'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.3);        
    plot(fold_in_eq(b_eq),fold_out_eq(b_eq),'Color',cmap(color_indices(i,2),:),'LineWidth',3);
    grid on
    box on    
    ylim([0 y_max])
    xlim([0 x_max])
    xlabel('log(c_{high}/c_{low})')
    ylabel('log(r_{high}/r_{low})')
    saveas(in_out_fold_eq,[FigPath type_names{i} '_fold-in_fold-out_eq.png'])
    saveas(in_out_fold_eq,[FigPath type_names{i} '_fold-in_fold-out_eq.pdf'])
    
    in_out_fold_neq = figure;
    hold on    
    scatter(fold_in_neq(plot_indices),fold_out_neq(plot_indices),mSize,'MarkerFaceColor',cmap(color_indices(i,1),:),'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.3);            
    scatter(fold_in_eq(plot_indices),fold_out_eq(plot_indices),mSize,'MarkerFaceColor',cmap(color_indices(i,2),:),'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.3);        
    pneq = plot(fold_in_neq(b_neq),fold_out_neq(b_neq),'Color',cmap(color_indices(i,1),:),'LineWidth',3);
    peq = plot(fold_in_eq(b_eq),fold_out_eq(b_eq),'Color',cmap(color_indices(i,2),:),'LineWidth',3);
    grid on
    box on  
    ylim([0 y_max])
    xlim([0 x_max])
    xlabel('log(c_{high}/c_{low})')
    ylabel('log(r_{high}/r_{low})')
    legend([pneq peq], 'non-equilibrium','equilibrium', 'Location','northwest','Fontsize',12) 
    saveas(in_out_fold_neq,[FigPath type_names{i} '_fold-in_fold-out_neq.png'])
    saveas(in_out_fold_neq,[FigPath type_names{i} '_fold-in_fold-out_neq.pdf'])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % next make basic sharpness fig    
    % extract vectors
    sharp_in_eq = reshape(metric_array_eq(:,3,:,:),[],1);
    sharp_out_eq = reshape(metric_array_eq(:,4,:,:),[],1);
    sharp_in_neq = reshape(metric_array_neq(:,3,:,:),[],1);
    sharp_out_neq = reshape(metric_array_neq(:,4,:,:),[],1);
    b_neq = convhull(sharp_in_neq,sharp_out_neq);
    b_eq = convhull(sharp_in_eq,sharp_out_eq);
    % draw indices to plot
    index_vec = 1:numel(sharp_out_neq);    
    plot_indices = randsample(index_vec,n_plot,false);
    % calculate plot boundaries
    y_max = max(sharp_out_neq(b_neq));
    x_max = max(sharp_in_neq(b_neq));
    % plot eq and non-eq scatters       a    
    in_out_sharp_eq = figure;
    hold on    
    scatter(sharp_in_eq(plot_indices),sharp_out_eq(plot_indices),mSize,'MarkerFaceColor',cmap(color_indices(i,2),:),'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.3);        
    plot(sharp_in_eq(b_eq),sharp_out_eq(b_eq),'Color',cmap(color_indices(i,2),:),'LineWidth',3);
    grid on
    box on    
    ylim([0 y_max])
    xlim([0 x_max])
    xlabel('\Delta c')
    ylabel('\Delta r')
    saveas(in_out_sharp_eq,[FigPath type_names{i} '_sharp-in_sharp-out_eq.png'])
    saveas(in_out_sharp_eq,[FigPath type_names{i} '_sharp-in_sharp-out_eq.pdf'])
    
    in_out_sharp_neq = figure;
    hold on    
    scatter(sharp_in_neq(plot_indices),sharp_out_neq(plot_indices),mSize,'MarkerFaceColor',cmap(color_indices(i,1),:),'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.3);            
    scatter(sharp_in_eq(plot_indices),sharp_out_eq(plot_indices),mSize,'MarkerFaceColor',cmap(color_indices(i,2),:),'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.3);        
    pneq = plot(sharp_in_neq(b_neq),sharp_out_neq(b_neq),'Color',cmap(color_indices(i,1),:),'LineWidth',3);
    peq = plot(sharp_in_eq(b_eq),sharp_out_eq(b_eq),'Color',cmap(color_indices(i,2),:),'LineWidth',3);
    grid on
    box on  
    ylim([0 y_max])
    xlim([0 x_max])
    xlabel('\Delta c')
    ylabel('\Delta r')
    legend([pneq peq], 'non-equilibrium','equilibrium', 'Location','northwest','Fontsize',12) 
    saveas(in_out_sharp_neq,[FigPath type_names{i} '_sharp-in_sharp-out_neq.png'])
    saveas(in_out_sharp_neq,[FigPath type_names{i} '_sharp-in_sharp-out_neq.pdf'])
%             
    metric_fig = figure;
    hold on                   
    % extract vectors
    mvec1_eq = reshape(metric_array_eq(:,metric_indices(1),:,:),[],1);
    mvec2_eq = reshape(metric_array_eq(:,metric_indices(2),:,:),[],1);
    mvec1_neq = reshape(metric_array_neq(:,metric_indices(1),:,:),[],1);
    mvec2_neq = reshape(metric_array_neq(:,metric_indices(2),:,:),[],1);
    
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

close all