clear 
close all

figPath = '../../param_sweeps/';
mkdir(figPath);

metric_names = {'fold-concentration','fold-production','fidelity','dynamicRange','dmdc','precision','decision_rate','error_rate'};
% generate all possible group sizes
base = 1:4;
edge_samp_cell = {base};
g_size = 1:4;
edge_options = 5:8;
for g = g_size
    c_array = nchoosek(edge_options,g);
    for c_vec = 1:size(c_array,1)
        edge_samp_cell = [edge_samp_cell{:} {[base c_array(c_vec,:)]}];
    end
end

index_array = [1,16; 
               1,16];
% call edge sampler for fidelity vs. absolute range
metric_indices = [3,8];
[sim_struct, fpath] = param_sweep(metric_indices,'edge_samp_cell', edge_samp_cell(unique(index_array(:))));

%%% make plots
type_names = {'activator', 'repressor'};
% plot metric dispersion

mSize = 15;
%%
close all
cmap = flipud(brewermap(9,'Set2'));
color_indices = [2,7,3];
for i = 1:numel(sim_struct)    
    metric_array = sim_struct(i).metric_array;
    sim_models = sim_struct.edge_samp_cell;
    type_models = edge_samp_cell(index_array(i,:));
    fig = figure;
    hold on    
    % iterate through different models
    p = [];  
    n_mdl = size(metric_array,5);
    iter = 1;
    for m = [1 2]%:n_mdl
        ind = n_mdl-m + 1;
        mdl = sim_models{ind};        
        mvec1 = reshape(metric_array(:,metric_indices(1),:,:,ind),[],1);
        mvec2 = reshape(metric_array(:,metric_indices(2),:,:,ind),[],1);
        
        % calculate boundaries        
        scatter(mvec1,-mvec2,mSize,'MarkerFaceColor',cmap(color_indices(iter),:),'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.3);        
        p = [p plot(0,0,'Color',cmap(color_indices(iter),:),'LineWidth',3)];  
        iter = iter + 1;
    end
    x_vec = linspace(0,.5);
    y_vec = 5 ./ x_vec;
    legend(p, 'neq full','eq','Location','best')
   
    grid on
    set(gca,'Fontsize',14)
    ylabel('decision rate')
    xlabel('dynamic range')
    saveas(fig,[figPath type_names{i} '_' metric_names{metric_indices(1)} '_' metric_names{metric_indices(2)} '.png'])
end

    