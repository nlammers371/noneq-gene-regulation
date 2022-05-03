% Plot results of sharpness parameter sweeps

% clear 
close all
addpath('../utilities/')

DataPath = '../../out/bivariate_parameter_sweeps_v3/';
FigPath = '../../fig/bivariate_parameter_sweeps_v3/';
mkdir(FigPath);

% get metrix names
[~,~, metric_names] = calculateMetricsGeneral([],[]);
flux_index = find(strcmp(metric_names,'Flux'));
tau_index = find(strcmp(metric_names,'CycleTime'));


%% (1) Plot Production Rate vs Sharpness
% set energy scale
atp = 30.5 * 1000 / 6.022e23; % in joules
kT = 300*1.38e-23;
n_plot = 1e3;

% set sim parameters
half_max_flag = 1;
metric_one_name = 'Sharpness';
metric_two_name = 'Precision';
metric_indices_pd_sharp = [find(strcmp(metric_names,metric_one_name)) find(strcmp(metric_names,metric_two_name))];

% generate load names
load_name_eq = ['param_sweep_results_' metric_one_name '_' ...
             metric_two_name '_eq1_hm' num2str(half_max_flag) '.mat'];
         
load_name_neq = ['param_sweep_results_' metric_one_name '_' ...
             metric_two_name '_eq0_hm' num2str(half_max_flag) '.mat'];    

         
% load data
load([DataPath load_name_eq]);
sim_results_sharp_prec_eq = simulation_results;
load([DataPath load_name_neq]);
sim_results_sharp_prec_neq = simulation_results;
clear simulation_results      
% make plots

% plot sharpness vs production rate (pos and neg sharpness)
mSize = 25;
close all

metric_array_eq = sim_results_sharp_prec_eq.metric_array;        
sharpness_vec_eq = metric_array_eq(:,metric_indices_pd_sharp(1));
precision_vec_eq = metric_array_eq(:,metric_indices_pd_sharp(2));

metric_array_neq = sim_results_sharp_prec_neq.metric_array;        
sharpness_vec_neq = metric_array_neq(:,metric_indices_pd_sharp(1));
precision_vec_neq = metric_array_neq(:,metric_indices_pd_sharp(2));

% calculate convex hull
% b_neq = convhull(mvec1_pd_sharp_eq,mvec2_pd );
plot_indices = randsample(1:length(sharpness_vec_eq),n_plot,false);
% plot eq and non-eq scatters           
pd_sharpness_full = figure;
cmap = brewermap([],'Set2');
colormap(cmap);
hold on    
scatter(sharpness_vec_neq,precision_vec_neq,mSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2, 'MarkerFaceColor',cmap(2,:));                
scatter(sharpness_vec_eq,precision_vec_eq,mSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2, 'MarkerFaceColor',cmap(3,:));                
grid on
box on    
% ylim([-.55 .55])
xlabel('sensitivity (dr/dc)')
ylabel('precision (1/\sigma_{int})')
set(gca,'FontSize',14)
% legend('non-equilibrium','equilibrium', 'Location','northeast') 
saveas(pd_sharpness_full,[FigPath metric_one_name '_' metric_two_name '_full.png'])
% saveas(pd_sharpness,[FigPath metric_names{metric_indices_pd(1)} '_' metric_names{metric_indices_pd(2)} '_full.pdf'])               


% plot eq only scatters           
pd_sharpness_act_eq = figure;
cmap = brewermap([],'Set2');
colormap(cmap);
hold on    
% scatter(abs(sharpness_vec_neq),precision_vec_neq,mSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2, 'MarkerFaceColor',cmap(2,:));                
scatter(abs(sharpness_vec_eq),precision_vec_eq,mSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2, 'MarkerFaceColor',cmap(3,:));                
grid on
box on    
xlabel('sensitivity (dr/dc)')
ylabel('precision (1/\sigma_{int})')
set(gca,'FontSize',14)
xlim([0 0.5])
% legend('non-equilibrium','equilibrium', 'Location','northeast') 
saveas(pd_sharpness_act_eq,[FigPath metric_one_name '_' metric_two_name '_act_eq.png'])


% plot eq and non-eq scatters           
pd_sharpness_act = figure;
cmap = brewermap([],'Set2');
colormap(cmap);
hold on    
scatter(abs(sharpness_vec_neq),precision_vec_neq,mSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2, 'MarkerFaceColor',cmap(2,:));                
scatter(abs(sharpness_vec_eq),precision_vec_eq,mSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2, 'MarkerFaceColor',cmap(3,:));                
grid on
box on    
xlabel('sensitivity (dr/dc)')
ylabel('precision (1/\sigma_{int})')
xlim([0 0.5])
set(gca,'FontSize',14)
% legend('non-equilibrium','equilibrium', 'Location','northeast') 
saveas(pd_sharpness_act,[FigPath metric_one_name '_' metric_two_name '_act.png'])

    
%% (2) Plot sharpness vs flux

% set sim parameters
half_max_flag = 1;
metric_one_name = 'Flux';
metric_two_name = 'Information';
metric_indices_info = [find(strcmp(metric_names,metric_one_name)) find(strcmp(metric_names,metric_two_name))];

% generate load names
load_name_neq = ['param_sweep_results_' metric_one_name '_' ...
             metric_two_name '_eq0_hm' num2str(half_max_flag) '.mat'];

load([DataPath load_name_neq]);


% estimate upper bound as a function of energy expenditure
metric_flux = simulation_results.metric_array;        
flux_vec = abs(metric_flux(:,metric_indices_info(1)))* kT / atp;
info_vec = abs(metric_flux(:,metric_indices_info(2)));

f_max = prctile(flux_vec,99.5);
flux_bins = linspace(-f_max,f_max);
f_window = 3*median(diff(flux_bins));

fs_vec_act = NaN(1,numel(flux_bins));

for  i = 1:numel(flux_bins)        
    f_ft = flux_vec >=flux_bins(i)-f_window & flux_vec <=flux_bins(i)+f_window;
    if any(f_ft)
        fs_vec_act(i) = nanmax(info_vec(f_ft));
    end
end

flux_info_fig = figure;
cmap2 = brewermap(9,'set2');
hold on
scatter(flux_vec,abs(info_vec),mSize,'MarkerfaceColor',cmap2(5,:),'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2);
p1 = plot(flux_bins,fs_vec_act,'--','Color','black','LineWidth',2);

xlim([0 f_max]);
% ylim([0.2 0.55])
% legend('network realizations','sharpness limit','Location','northwest')

xlabel('energy dissipation per cycle (ATP units)')
ylabel('idealized information (dr/dc  \sigma_{int}^{-1})')
grid on
set(gca,'FontSize',14)

saveas(flux_info_fig,[FigPath metric_one_name '_' metric_two_name '.png'])

%% (3) Look at interplay between precision, sharpness, and cycle time
kappa_vec_eq = 1./metric_array_eq(:,tau_index);
kappa_vec_neq = 1./metric_array_neq(:,tau_index);
kappa_vec_flux = 1./metric_flux (:,tau_index);
close all

% plot eq and non-eq scatters           
pd_sharpness_act_tau = figure;
cmap1 = flipud(brewermap([],'Spectral'));
colormap(cmap1);
hold on    
scatter(abs(sharpness_vec_neq),precision_vec_neq,mSize,kappa_vec_neq,'filled','MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2);%, 'MarkerFaceColor',cmap(2,:));                
scatter(abs(sharpness_vec_eq),precision_vec_eq,mSize,kappa_vec_eq,'filled','MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2);%, 'MarkerFaceColor',cmap(3,:));                
grid on
box on    
xlabel('sensitivity (dr/dc)')
ylabel('precision (1/\sigma_{int})')
h = colorbar;
ylabel(h,'cycle rate (s^{-1})')
set(gca,'FontSize',14)
saveas(pd_sharpness_act_tau,[FigPath 'Sharpness_Precision_tau.png'])

%%
pd_sharpness_tau_eq = figure;
cmap = brewermap([],'Set2');
colormap(cmap);
hold on    
% scatter(kappa_vec_neq,abs(sharpness_vec_neq),mSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2, 'MarkerFaceColor',cmap(2,:));                
scatter(kappa_vec_eq(plot_indices),abs(sharpness_vec_eq(plot_indices)),mSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2, 'MarkerFaceColor',cmap(3,:));                
grid on
box on    
ylabel('sensitivity (dr/dc)')
xlabel('cycle rate (s^{-1})')
set(gca,'FontSize',14)
set(gca,'xScale','log')
% legend('non-equilibrium','equilibrium', 'Location','northeast') 
saveas(pd_sharpness_tau_eq,[FigPath 'Sharpness_CycleTime_eq.png'])
saveas(pd_sharpness_tau_eq,[FigPath 'Sharpness_CycleTime_eq.pdf'])

pd_sharpness_tau = figure;
cmap = brewermap([],'Set2');
colormap(cmap);
hold on    
scatter(kappa_vec_neq(plot_indices),abs(sharpness_vec_neq(plot_indices)),mSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2, 'MarkerFaceColor',cmap(2,:));                
scatter(kappa_vec_eq(plot_indices),abs(sharpness_vec_eq(plot_indices)),mSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2, 'MarkerFaceColor',cmap(3,:));                
grid on
box on    
ylabel('sensitivity (dr/dc)')
xlabel('cycle rate (s^{-1})')
set(gca,'FontSize',14)
set(gca,'xScale','log')
% legend('non-equilibrium','equilibrium', 'Location','northeast') 
saveas(pd_sharpness_tau,[FigPath 'Sharpness_CycleTime.png'])
saveas(pd_sharpness_tau,[FigPath 'Sharpness_CycleTime.pdf'])

pd_sharpness_tau_neq = figure;
cmap = brewermap([],'Set2');
colormap(cmap);
hold on    
scatter(kappa_vec_neq(plot_indices),abs(sharpness_vec_neq(plot_indices)),mSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2, 'MarkerFaceColor',cmap(2,:));                
grid on
box on    
ylabel('sensitivity (dr/dc)')
xlabel('cycle rate (s^{-1})')
set(gca,'FontSize',14)
ylim([0 0.55])
set(gca,'xScale','log')
% legend('non-equilibrium','equilibrium', 'Location','northeast') 
saveas(pd_sharpness_tau_neq,[FigPath 'Sharpness_CycleTime_neq.png'])
saveas(pd_sharpness_tau_neq,[FigPath 'Sharpness_CycleTime_neq.pdf'])

%%
pd_precision_tau_eq = figure;
cmap = brewermap([],'Set2');
colormap(cmap);
hold on    
% scatter(kappa_vec_neq,precision_vec_neq,mSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2, 'MarkerFaceColor',cmap(2,:));                
scatter(kappa_vec_eq,precision_vec_eq,mSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2, 'MarkerFaceColor',cmap(3,:));                
grid on
box on    
ylabel('precision (1/\sigma_{int})')
xlabel('cycle rate (s^{-1})')
set(gca,'FontSize',14)
set(gca,'xScale','log')
set(gca,'yScale','log')
% legend('non-equilibrium','equilibrium', 'Location','northeast') 
saveas(pd_precision_tau_eq,[FigPath 'Precision_CycleTime_eq.png'])


pd_precision_tau = figure;
cmap = brewermap([],'Set2');
colormap(cmap);
hold on    
scatter(kappa_vec_neq,precision_vec_neq,mSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2, 'MarkerFaceColor',cmap(2,:));                
scatter(kappa_vec_eq,precision_vec_eq,mSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2, 'MarkerFaceColor',cmap(3,:));                
grid on
box on    
ylabel('precision (1/\sigma_{int})')
xlabel('cycle rate (s^{-1})')
set(gca,'FontSize',14)
set(gca,'xScale','log')
set(gca,'yScale','log')
% legend('non-equilibrium','equilibrium', 'Location','northeast') 
saveas(pd_precision_tau,[FigPath 'Precision_CycleTime.png'])

