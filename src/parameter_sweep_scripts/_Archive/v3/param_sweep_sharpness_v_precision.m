clear 
close all

figPath = '../../fig/param_sweep/';
mkdir(figPath);

metric_names = {'fold-concentration','fold-production','fidelity','dynamicRange','dmdc','precision'};
% call edge sampler for fidelity vs. absolute range
metric_indices = [5,6];
[sim_struct, fpath] = param_sweep(metric_indices);

%%% make plots
type_names = {'activator', 'repressor'};
% plot metric dispersion

mSize = 15;
%%
close all
cmap = flipud(brewermap(9,'Set2'));
color_indices = [2,3];
plot_indices = [4 1];
for i = 1:numel(sim_struct)    
    metric_array = sim_struct(i).metric_array;
    sim_models = sim_struct.edge_samp_cell;    
    fig = figure;
    hold on    
    % iterate through different models
    p = [];  
    n_mdl = size(metric_array,5);
    iter = 1;
    for m = plot_indices        
        % extract vectors
        mvec1 = reshape(metric_array(:,metric_indices(1),:,:,m),[],1);
        mvec2 = reshape(metric_array(:,metric_indices(2),:,:,m),[],1);        
        % remove low sharpness
        ft = mvec1 > .02;
        mvec2 = mvec2(ft);
        mvec1 = mvec1(ft);                        
        mvec2 = exp(mvec2);        
        % calculate boundaries        
        scatter(mvec1,mvec2,mSize,'MarkerFaceColor',cmap(color_indices(iter),:),'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.3);        
        p = [p plot(0,0,'Color',cmap(color_indices(iter),:),'LineWidth',3)];  
        iter = iter + 1;
    end
    xlim([.05 .5])    
    ylim([0 10])
  
    legend(p, 'non-equilibrium', 'equilibrium','Location','best')
   
    grid on
    set(gca,'Fontsize',14)
    ylabel('precision')
    xlabel('sharpness')
    saveas(fig,[figPath type_names{i} '_' metric_names{metric_indices(1)} '_' metric_names{metric_indices(2)} '.png'])
    
    % plot using color to convey info rate
    fig = figure;
    colormap(flipud(brewermap(128,'Spectral')));
    hold on        
    n_mdl = size(metric_array,5);
    iter = 1;
    for m = plot_indices        
        % extract vectors
        mvec1 = reshape(metric_array(:,metric_indices(1),:,:,m),[],1);
        mvec2 = reshape(metric_array(:,metric_indices(2),:,:,m),[],1);        
        % remove low sharpness
        ft = mvec1 > .02;
        mvec2 = mvec2(ft);
        mvec1 = mvec1(ft);                        
        mvec2 = exp(mvec2);        
        % calculate boundaries        
        scatter(mvec1,mvec2,mSize,mvec1.*mvec2,'filled','MarkerEdgeAlpha',0,'MarkerFaceAlpha',.3);                
        iter = iter + 1;
    end
    xlim([.05 .5])    
    ylim([0 10])
    h = colorbar;    
    grid on    
    ylabel('precision')
    ylabel(h,'information rate')
    xlabel('sharpness')
    set(gca,'Fontsize',14)
    saveas(fig,[figPath type_names{i} '_' metric_names{metric_indices(1)} '_' metric_names{metric_indices(2)} 'info_hm.png'])
end

    