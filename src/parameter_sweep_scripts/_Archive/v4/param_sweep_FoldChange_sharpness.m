% script to generate data and figures for fidelity v sharpness parameter
% space exploration

clear 
close all
addpath('../utilities/')
FigPath = '../../fig/bivariate_param_sweeps/';
mkdir(FigPath);
[~, metric_names] = calculateMetricsV2([],[],[]);

% call edge sampler for fidelity vs. absolute range
metric_indices = [3,4];
quadrant_ref = [1,1];
quadrant_tol = .97;
% [sim_struct_eq, fpath_eq] = param_sweep(metric_indices,1);
% [sim_struct_neq, fpath_neq] = param_sweep(metric_indices,0);

metric_indices_noise = [4,6];
% [sim_struct_ns_eq, ~] = param_sweep(metric_indices_noise,1);%,'n_seeds',10);
[sim_struct_ns_neq, ~] = param_sweep(metric_indices_noise,0,'quadrant_tol',quadrant_tol);

%% make plots
type_names = {'burst frequency', 'burst duration'};
% plot metric dispersion
mSize = 25;
close all
cmap1 = brewermap(9,'Set2');
cmap2 = brewermap(128,'Spectral');
% n_plot = 2000;
cmap = brewermap(128,'Spectral');
for i = 2%1:numel(sim_struct_ns_eq)        
    close all
    if i == 1
        cn = cmap2(115,:);
        ce = cmap2(90,:);
    else
        cn = cmap1(3,:);
        ce = cmap1(2,:);
    end
    %%% fold change vs. sharpness
    metric_array_eq = sim_struct_eq(i).metric_array;
    metric_array_neq = sim_struct_neq(i).metric_array;            
    
    metric_fig_eq = figure;
    hold on                   
    % extract vectors
    mvec1_eq = reshape(metric_array_eq(:,metric_indices(1),:,:),[],1);
    mvec2_eq = reshape(metric_array_eq(:,metric_indices(2),:,:),[],1);
    mvec1_neq = reshape(metric_array_neq(:,metric_indices(1),:,:),[],1);
    mvec2_neq = reshape(metric_array_neq(:,metric_indices(2),:,:),[],1);
    
    % plot eq and non-eq scatters       
%     scatter(mvec1_neq,mvec2_neq,mSize,'MarkerFaceColor',cn,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2);        
    scatter(mvec1_eq,mvec2_eq,mSize,'MarkerFaceColor',ce,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2);        

    grid on
    box on    
    xlabel('normalized fold-change')
    ylabel('sharpness')
    set(gca,'Fontsize',14)
    axis([0 2 0 .5])
%     legend('non-equilibrium','equilibrium', 'Location','northeast','Fontsize',12) 
    saveas(metric_fig_eq,[FigPath type_names{i} '_' metric_names{metric_indices(1)} '_' metric_names{metric_indices(2)} '_eq_only.png'])
    
    metric_fig_neq = figure;
    hold on                   
    
    % plot eq and non-eq scatters       
    scatter(mvec1_neq,mvec2_neq,mSize,'MarkerFaceColor',cn,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2);        
    scatter(mvec1_eq,mvec2_eq,mSize,'MarkerFaceColor',ce,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2);        

    grid on
    box on    
    xlabel('normalized fold-change')
    ylabel('sharpness')
    set(gca,'Fontsize',14)
%     legend('non-equilibrium','equilibrium', 'Location','northeast','Fontsize',12) 
    saveas(metric_fig_neq,[FigPath type_names{i} '_' metric_names{metric_indices(1)} '_' metric_names{metric_indices(2)} 'neq.png'])


%%% noise vs. sharpness
%     metric_array_ns_eq = sim_struct_ns_eq(i).metric_array;
%     metric_array_ns_neq = sim_struct_ns_neq(i).metric_array;            
%     
%     metric_fig_eq = figure;
%     hold on                   
%     % extract vectors
%     mvec1_ns_eq = reshape(metric_array_ns_eq(:,metric_indices_noise(1),:,:),[],1);
%     mvec2_ns_eq = reshape(metric_array_ns_eq(:,metric_indices_noise(2),:,:),[],1);
%     mvec1_ns_neq = reshape(metric_array_ns_neq(:,metric_indices_noise(1),:,:),[],1);
%     mvec2_ns_neq = reshape(metric_array_ns_neq(:,metric_indices_noise(2),:,:),[],1);
%     
%     % plot eq and non-eq scatters       
% %     scatter(mvec1_ns_neq,mvec2_ns_neq,mSize,'MarkerFaceColor',cn,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2);        
%     scatter(mvec1_ns_eq,mvec2_ns_eq,mSize,'MarkerFaceColor',ce,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2);        
% 
%     grid on
%     axis([.05,.5 0 15])
%     box on    
%     ylabel('precision (1/\sigma_r)')
%     xlabel('sharpness (dr/dc)')
%     set(gca,'Fontsize',14)
% %     legend('non-equilibrium','equilibrium', 'Location','northeast','Fontsize',12) 
%     saveas(metric_fig_eq,[FigPath type_names{i} '_' metric_names{metric_indices_noise(1)} '_' metric_names{metric_indices_noise(2)} '_eq_only.png'])
%     
%     metric_fig_neq = figure;
%     hold on                           
%     % plot eq and non-eq scatters       
%     scatter(mvec1_ns_neq,mvec2_ns_neq,mSize,'MarkerFaceColor',cn,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2);        
%     scatter(mvec1_ns_eq,mvec2_ns_eq,mSize,'MarkerFaceColor',ce,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2);        
% 
%     grid on
%     axis([.05,.5 0 15])
%     box on    
%     ylabel('precision (1/\sigma_r)')
%     xlabel('sharpness (dr/dc)')
%     set(gca,'Fontsize',14)
% %     legend('non-equilibrium','equilibrium', 'Location','northeast','Fontsize',12) 
%     saveas(metric_fig_neq,[FigPath type_names{i} '_' metric_names{metric_indices_noise(1)} '_' metric_names{metric_indices_noise(2)} '_neq.png'])
end

    