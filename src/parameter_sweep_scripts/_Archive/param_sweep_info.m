% script to generate data and figures for fidelity v sharpness parameter
% space exploration
clear 
close all
addpath('../utilities/')
FigPath = '../../fig/bivariate_param_sweeps_v2/';
mkdir(FigPath);
[~, metric_names] = calculateMetricsV3([],[]);

% call edge sampler for production rate vs. sharpness
metric_indices_var = [1,4];
metric_indices_flux = [3,5];
[sim_struct_flux, ~] = param_sweep_v2(metric_indices_flux,'n_seeds',5);
[sim_struct_var, ~] = param_sweep_v2(metric_indices_var,'n_seeds',5);


% set energy scale
atp = 30.5 * 1000 /6.022e23; % in joules
kT = 1.38e-23*300;

%% make plots
% plot metric dispersion
type_names = {'activator','repressor'};
mSize = 25;

close all
cmap = flipud(brewermap([],'Spectral'));
n_plot = 5e3;

for i = 1:numel(sim_struct_var)        
    %%%%%%%%%% production rate vs. sharpness
    metric_array_pd_neq = sim_struct_var(i).metric_array;        
    mvec1_pd_neq = metric_array_pd_neq(:,metric_indices_var(1));
    mvec2_pd_neq = metric_array_pd_neq(:,metric_indices_var(2));
    flux_pd_neq = metric_array_pd_neq(:,3) * kT / atp;
    b_neq = convhull(mvec1_pd_neq,mvec2_pd_neq );
   
    % draw indices to plot
    index_vec = 1:numel(mvec1_pd_neq);      
    plot_indices = randsample(index_vec,n_plot,false);
    % plot eq and non-eq scatters           
    pd_sharpness = figure;
    colormap(cmap);
    hold on    
    scatter(mvec1_pd_neq(plot_indices),mvec2_pd_neq(plot_indices),mSize,flux_pd_neq(plot_indices),'filled','MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2);                
    grid on
    box on    
    h = colorbar;
    xlabel('normalized mRNA production rate')
    ylabel('sharpness')
    ylabel(h,'energy dissipation (ATP units)')
    set(gca,'FontSize',14)
%     legend([pneq peq], 'non-equilibrium','equilibrium', 'Location','northeast','Fontsize',12) 
    saveas(pd_sharpness,[FigPath type_names{i} '_' metric_names{metric_indices_var(1)} '_' metric_names{metric_indices_var(2)} '.png'])
    saveas(pd_sharpness,[FigPath type_names{i} '_' metric_names{metric_indices_var(1)} '_' metric_names{metric_indices_var(2)} '.pdf'])               
end

    
%% Plot sharpness vs flux

metric_act = sim_struct_flux(1).metric_array;        
flux_vec_act = metric_act(:,3)* kT / atp;
info_vec_act = metric_act(:,5);

metric_rep = sim_struct_flux(2).metric_array;        
flux_vec_rep = metric_rep(:,3)* kT / atp;
info_vec_rep = metric_rep(:,5);

f_max = max(abs(vertcat(flux_vec_rep, flux_vec_act)));
flux_bins = linspace(0,f_max);
f_window = 1*median(diff(flux_bins));

fs_vec_act = NaN(1,numel(flux_bins));
fs_vec_rep = NaN(1,numel(flux_bins));

for  i = 1:numel(flux_bins)        
    f_ft_act = flux_vec_act >=flux_bins(i)-f_window & flux_vec_act <=flux_bins(i)+f_window;
    fs_vec_act(i) = nanmax(abs(info_vec_act(f_ft_act)));
    f_ft_rep = flux_vec_rep >=flux_bins(i)-f_window & flux_vec_rep <=flux_bins(i)+f_window;
    fs_vec_rep(i) = nanmax(abs(info_vec_rep(f_ft_rep)));
end

plot_indices = randsample(1:numel(info_vec_act),n_plot,false);

close all
flux_info_fig = figure;
cmap2 = brewermap(9,'set2');
hold on
scatter(flux_vec_act(plot_indices),info_vec_act(plot_indices),mSize,'MarkerfaceColor',cmap2(2,:),'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2);
scatter(flux_vec_rep(plot_indices),info_vec_rep(plot_indices),mSize,'MarkerfaceColor',cmap2(3,:),'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2);
p2 = plot(flux_bins,imgaussfilt(fs_vec_rep,2),'-','Color',cmap2(3,:),'LineWidth',2);
p1 = plot(flux_bins,imgaussfilt(fs_vec_act,2),'--','Color',cmap2(2,:),'LineWidth',2);

xlim([0 f_max]);
ylim([1e-5 10*(max([fs_vec_rep fs_vec_act]))])
legend([p1 p2],'activator','repressor','Location','southwest')
% set(gca,'YScale','log')

xlabel('energy dissipation (ATP units)')
ylabel('information rate')

set(gca,'YScale','log')
set(gca,'FontSize',14)
grid on

% 
saveas(flux_info_fig,[FigPath 'info_vs_flux_plot.png'])
saveas(flux_info_fig,[FigPath 'info_vs_flux_plot.pdf'])

%%
act_rates = sim_struct_flux(1).rate_array;
[~,sort_indices] = sort(info_vec_act);
[max_info,max_index] = max(info_vec_act);