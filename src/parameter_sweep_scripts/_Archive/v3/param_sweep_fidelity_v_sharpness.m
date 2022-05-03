% script to generate data and figures for fidelity v sharpness parameter
% space exploration

clear 
close all

figPath = '../../fig/param_sweep/';
mkdir(figPath);
metric_names = {'fold-concentration','fold-production','fidelity','dynamicRange','dmdc','precision'};

% call edge sampler for fidelity vs. absolute range
metric_indices = [3,4];
[sim_struct, fpath] = param_sweep(metric_indices);

%%% make plots
type_names = {'activator', 'repressor'};
% plot metric dispersion
mSize = 15;
close all
cmap = flipud(brewermap(9,'Set2'));
color_indices = [2,7,8,5,3];
%%
for i = 1:numel(sim_struct)    
    activator_flag = i==1;
    metric_array = sim_struct(i).metric_array;
    sim_models = sim_struct.edge_samp_cell;    
    fig = figure;
    hold on    
    % iterate through different models
    p = [];  
    n_mdl = size(metric_array,5);
    iter = 1;
    for m = 1:n_mdl
        % jenk
        if m == n_mdl-1
            iter = iter + 1;
            continue
        end
        ind = n_mdl-m + 1;        
        % extract vectors
        mvec1 = reshape(metric_array(:,metric_indices(1),:,:,ind),[],1);
        mvec2 = reshape(metric_array(:,metric_indices(2),:,:,ind),[],1);
                
        % calculate boundaries        
        scatter(mvec1,mvec2,mSize,'MarkerFaceColor',cmap(color_indices(iter),:),'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.3);        
        p = [p plot(0,0,'Color',cmap(color_indices(iter),:),'LineWidth',3)];  
        iter = iter + 1;
    end      
    grid on
    set(gca,'Fontsize',14)
    xlabel('fold-change')
    ylabel('sharpness')
    legend(p, 'neq full','neq (non-specific)', 'neq (specific)', 'equilibrum', ...
        'Location','best','Fontsize',12) 
    saveas(fig,[figPath type_names{i} '_' metric_names{metric_indices(1)} '_' metric_names{metric_indices(2)} '.png'])
end

    