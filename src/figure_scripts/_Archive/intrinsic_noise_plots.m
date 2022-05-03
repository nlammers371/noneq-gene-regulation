% Plot results of sharpness vs. intrinsic noise parameter sweeps

clear 
close all
addpath(genpath('../utilities/'))

DataPath = '../../out/bivariate_parameter_sweeps/';
FigPath = '../../fig/intrinsic_noise_plots/';
mkdir(FigPath);

% get metric names
[~,metric_names] = calculateMetrics_v4([]);
flux_index = find(strcmp(metric_names,'Flux'));


% set energy scale
atp = 50 * 1000 / 6.022e23; % in joules
kbT = 300*1.38e-23;

% specify number of dots to plot
n_plot = 4e3;

% set plot colors
pboc = [228 221 210]/256;

% set sim parameters
metric_indices_int = [1,4];
metric_one_name = metric_names{metric_indices_int(1)};
metric_two_name = metric_names{metric_indices_int(2)};

% generate load names
load_name_eq = ['param_sweep_results_' metric_names{metric_indices_int(1)} '_' ...
             metric_names{metric_indices_int(2)} '_eq1.mat'];
         
load_name_neq = ['param_sweep_results_' metric_names{metric_indices_int(1)} '_' ...
             metric_names{metric_indices_int(2)} '_eq0.mat'];    

         
% load data
load([DataPath load_name_eq]);
load([DataPath load_name_neq]);    
% make plots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot sensitivity vs intrinsic noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mSize = 25;
close all

metric_array_eq = sim_struct_int_eq.metric_array;        
mvec1_int_sharp_eq = reshape(metric_array_eq(:,metric_indices_int(1),:),1,[]);
mvec2_int_sharp_eq = reshape(metric_array_eq(:,metric_indices_int(2),:),1,[]);

metric_array_neq = sim_struct_int.metric_array;        
mvec1_int_sharp_neq = reshape(metric_array_neq(:,metric_indices_int(1),:),1,[]);
mvec2_int_sharp_neq = reshape(metric_array_neq(:,metric_indices_int(2),:),1,[]);

% esitmate boundary and resample points to emphasize those near limits for
% each condition 
% b_neq = boundary([mvec1_int_sharp_neq',log(mvec2_int_sharp_neq')],1);
% bv1_neq = mvec1_int_sharp_neq(b_neq);
% bv2_neq = mvec2_int_sharp_neq(b_neq);
% dist_vec = min(sqrt((mvec1_int_sharp_neq' - bv1_neq).^2+(mvec2_int_sharp_neq' - bv2_neq).^2),[],2);

plot_indices_neq = randsample(find(mvec1_int_sharp_neq>=0),n_plot,false);
plot_indices_eq = randsample(find(mvec1_int_sharp_eq>=0),n_plot,false);
%%
% figure with full model behavior (activator and repressor) 
sensitivity_int = figure;
cmap = brewermap(9,'Set2');
colormap(cmap);
hold on    
scatter(mvec1_int_sharp_neq(plot_indices_neq),mvec2_int_sharp_neq(plot_indices_neq),mSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2, 'MarkerFaceColor',cmap(2,:));                
scatter(mvec1_int_sharp_eq(plot_indices_eq),mvec2_int_sharp_eq(plot_indices_eq),mSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2, 'MarkerFaceColor',cmap(3,:));                
grid on
box on    
% ylim([1e-1 1e4])
xlabel('sensitivity (dr/dc)')
ylabel('precision (\sigma_{int}^{-1})')
set(gca,'FontSize',14)
% set(gca,'yScale','log','xtick',-0.5:0.25:0.5)
set(gca,'Color',pboc)
% legend('non-equilibrium','equilibrium', 'Location','northeast') 

sensitivity_int.InvertHardcopy = 'off';
saveas(sensitivity_int,[FigPath metric_one_name '_' metric_two_name '.png'])
saveas(sensitivity_int,[FigPath metric_one_name '_' metric_two_name '.pdf'])



%% Information vs flux
% set sim parameters
close all
metric_one_name = 'Flux';
metric_two_name = 'Information';
metric_indices_pd_sharp = [find(strcmp(metric_names,metric_one_name)) find(strcmp(metric_names,metric_two_name))];

% generate load names
load_name_neq = ['param_sweep_results_' metric_one_name '_' ...
             metric_two_name '_eq0.mat'];
% load
load([DataPath load_name_neq]);

% extract vectors
metric_flux = sim_struct_flux.metric_array;        
flux_vec = reshape(metric_flux(:,flux_index,:),1,[])* kbT / atp;
info_vec = reshape(metric_flux(:,5,:),1,[]);
f_max = prctile(flux_vec,99.9);

% generate plot indices
plot_indices_flux = randsample(find(info_vec>0&flux_vec>0),length(find(info_vec>0&flux_vec>0)),false);
% calculate boundary
% b_neq = boundary([flux_vec(plot_indices_flux)',info_vec(plot_indices_flux)'],1);


flux_info_fig = figure;
cmap2 = brewermap(9,'set2');
hold on
scatter(flux_vec(plot_indices_flux),info_vec(plot_indices_flux),mSize,'MarkerfaceColor',cmap2(5,:),'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.4);
% p1 = plot(flux_bins,fs_vec_act,'-','Color','black','LineWidth',2);

xlim([0 1.2]);
ylim([0 1.1]);
set(gca,'Color',pboc)
xlabel('energy dissipation per cycle (ATP units)')
ylabel('information per cycle')
% grid on
set(gca,'FontSize',14)
box on
flux_info_fig.InvertHardcopy = 'off';
saveas(flux_info_fig,[FigPath 'info_vs_flux_plot.png'])
saveas(flux_info_fig,[FigPath 'info_vs_flux_plot.pdf'])
