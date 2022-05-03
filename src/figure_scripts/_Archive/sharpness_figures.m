% Plot results of sharpness parameter sweeps
clear 
close all
addpath(genpath('../utilities/'))

% %%%%%%%%%%%%%%%%  set relative read and write paths %%%%%%%%%%%%%%%%%%%%
DataPath = '../../out/bivariate_parameter_sweeps/';
FigPath = '../../fig/sharpness_plots/';
mkdir(FigPath);

% %%%%%%%%%%%%%%%%  get metric names %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,metric_names] = calculateMetrics_v4([]);
flux_index = find(strcmp(metric_names,'Flux'));
rate_index = find(strcmp(metric_names,'Production Rate'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
frac_index = find(strcmp(metric_names,'2StateFrac'));
% %%%%%%%%%%%%%%%%  Set plot parameters and constants %%%%%%%%%%%%%%%%%%%%

atp = 50 * 1000 / 6.022e23; % in joules
kbT = 300*1.38e-23;


n_plot = 3e3; % number of points to plot
markerAlpha = 0.5; % marker transparency
markerSize = 25; % marker size

%% %%%%%%%%%%%%%%%%  Plot Production Rate vs Sharpness %%%%%%%%%%%%%%%%%%%%

% set sim parameters
s_v_r_indices = [sharpness_index rate_index];
metric_one_name = metric_names{s_v_r_indices(1)} ;
metric_two_name = metric_names{s_v_r_indices(2)} ;

% generate load names
load_name_eq = ['param_sweep_results_' metric_one_name '_' metric_two_name ...
             '_eq1.mat'];
         
load_name_neq = ['param_sweep_results_' metric_one_name '_' metric_two_name ...
              '_eq0.mat'];    
         
% load data
load([DataPath load_name_eq]);
load([DataPath load_name_neq]);    

% make plots
close all

% extract parameter vectors
metric_array_eq = sim_struct_pd_eq.metric_array;        
mvec1_rate_v_sharp_eq = reshape(metric_array_eq(:,s_v_r_indices(1),:),1,[]);
mvec2_rate_v_sharp_eq = reshape(metric_array_eq(:,s_v_r_indices(2),:),1,[]);

metric_array_neq = sim_struct_pd.metric_array;        
mvec1_rate_v_sharp_neq = reshape(metric_array_neq(:,s_v_r_indices(1),:),1,[]);
mvec2_rate_v_sharp_neq = reshape(metric_array_neq(:,s_v_r_indices(2),:),1,[]);

plot_indices_eq = randsample(find(mvec1_rate_v_sharp_eq>=0),n_plot,false);
plot_indices_neq = randsample(find(mvec1_rate_v_sharp_neq>=0),n_plot,false);

%% plot eq and non-eq scatters           
pd_sharpness_basic = figure;

cmap = brewermap(9,'Set2');

hold on    
scatter(mvec2_rate_v_sharp_neq(plot_indices_neq),mvec1_rate_v_sharp_neq(plot_indices_neq),...
  markerSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(2,:));                
scatter(mvec2_rate_v_sharp_eq(plot_indices_eq),mvec1_rate_v_sharp_eq(plot_indices_eq),...
  markerSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(3,:));                

box on    

StandardFigurePBoC([],gca)

ylim([0 .52])
xlabel('production rate (p_{on})')
ylabel('sharpness (dr/ dc)')
set(gca,'FontSize',14)

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
set(gca,'YTick',0:.1:.5)
set(gca,'XTick',0:.25:1)
pd_sharpness_basic.InvertHardcopy = 'off';


saveas(pd_sharpness_basic,[FigPath metric_one_name '_' metric_two_name '.png'])
saveas(pd_sharpness_basic,[FigPath metric_one_name '_' metric_two_name '.pdf'])

% %%%%%%%%%%%%%%%%  Add simple 2 state predictions %%%%%%%%%%%%%%%%%%%%%%%  
r_vec = linspace(0,1);
pd_neq_lim = 2*(1-r_vec).*r_vec;
pd_eq_lim = (1-r_vec).*r_vec;

pd_sharpness_pd = figure;
colormap(cmap);
hold on    
scatter(mvec2_rate_v_sharp_neq(plot_indices_neq),mvec1_rate_v_sharp_neq(plot_indices_neq),...
  markerSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(2,:));                
scatter(mvec2_rate_v_sharp_eq(plot_indices_eq),mvec1_rate_v_sharp_eq(plot_indices_eq),...
  markerSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(3,:)); 

plot(r_vec,pd_neq_lim,'-','Color','black','LineWidth',1.5)
plot(r_vec,pd_eq_lim,'--','Color','black','LineWidth',1.5)


StandardFigurePBoC([],gca)

% grid on
box on    
ylim([0 .52])

xlabel('production rate (p_{on})')
ylabel('sharpness (dr/dc)')
set(gca,'FontSize',14)

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
set(gca,'YTick',0:.1:.5)
set(gca,'XTick',0:.25:1)
pd_sharpness_pd.InvertHardcopy = 'off';

% legend('non-equilibrium','equilibrium', 'Location','northeast') 
saveas(pd_sharpness_pd,[FigPath metric_one_name '_' metric_two_name '_theory.png'])
saveas(pd_sharpness_pd,[FigPath metric_one_name '_' metric_two_name '_theory.pdf'])


%% %%%%%%%%%% Plot on and off rate c dependencies %%%%%%%%%%%%%%%%%%%%%%%%%
close all

% get variable indices
kon_index = find(strcmp(metric_names,'KonEff'));
koff_index = find(strcmp(metric_names,'KoffEff'));
dkon_index = find(strcmp(metric_names,'dKondC'));
dkoff_index = find(strcmp(metric_names,'dKoffdC'));
dkdc_indices = [dkon_index dkoff_index];
% name of data sets
save_name_dkdc = ['param_sweep_results_' metric_names{dkdc_indices(1)} '_' ...
    metric_names{dkdc_indices(2)}];

% load
eq_flags = [0 1];
dkdc_struct = struct;
for e = 1:length(eq_flags)
  load([DataPath save_name_dkdc '_eq' num2str(eq_flags(e)) '.mat'],'sim_struct')
  dkdc_struct(e).sim_struct = sim_struct;
end

% off and on rate c dependencies
koffNeq = reshape(dkdc_struct(1).sim_struct.metric_array(:,koff_index,:),1,[]);
konNeq = reshape(dkdc_struct(1).sim_struct.metric_array(:,kon_index,:),1,[]);
dKoffdCNeq = reshape(dkdc_struct(1).sim_struct.metric_array(:,dkoff_index,:),1,[]);%./koffNeq;
dKondCNeq = reshape(dkdc_struct(1).sim_struct.metric_array(:,dkon_index,:),1,[]);%./konNeq;

pdRateEq = reshape(dkdc_struct(2).sim_struct.metric_array(:,2,:),1,[]);
sharpnessEq = reshape(dkdc_struct(2).sim_struct.metric_array(:,1,:),1,[])./(pdRateEq.*(1-pdRateEq));
fracEq = reshape(dkdc_struct(2).sim_struct.metric_array(:,frac_index,:),1,[]);

koffEq = reshape(dkdc_struct(2).sim_struct.metric_array(:,koff_index,:),1,[]);
konEq = reshape(dkdc_struct(2).sim_struct.metric_array(:,kon_index,:),1,[]);
dKoffdCEq = reshape(dkdc_struct(2).sim_struct.metric_array(:,dkoff_index,:),1,[]);%./koffEq;
dKondCEq = reshape(dkdc_struct(2).sim_struct.metric_array(:,dkon_index,:),1,[]);%./konEq;

plot_indices_eq = randsample(find(sharpnessEq>=0.0),min([n_plot length(find(sharpnessEq>=0.0))]),false);

fluxNeq = abs(reshape(dkdc_struct(1).sim_struct.metric_array(:,3,:),1,[]));
pdRateNeq = reshape(dkdc_struct(1).sim_struct.metric_array(:,2,:),1,[]);
sharpnessNeq = reshape(dkdc_struct(1).sim_struct.metric_array(:,1,:),1,[])./(pdRateNeq.*(1-pdRateNeq));

plot_indices_neq = randsample(find(sharpnessNeq>=0.0),min([n_plot length(find(sharpnessNeq>=0.0))]),false);


% make figure
flux_rate_fig = figure;
cmap3 = flipud(brewermap([],'Spectral'));
hold on
colormap(cmap3);

scatter(dKondCNeq(plot_indices_neq),dKoffdCNeq(plot_indices_neq),markerSize,sharpnessNeq(plot_indices_neq),'filled',...
        'MarkerEdgeAlpha',.1,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha);

scatter(dKondCEq(plot_indices_eq),dKoffdCEq(plot_indices_eq),markerSize,sharpnessEq(plot_indices_eq),'filled',...
        'MarkerEdgeAlpha',.1,'MarkerFaceAlpha',markerAlpha,'MarkerEdgeColor','k');


h = colorbar;
ylabel(h,'normalized sharpness')
caxis([0 2])
xlim([0 2.1]);
ylim([-1 0]);
% ylim([-1 0])
xlabel('d k_{on} / dc')
ylabel('d k_{off} / dc')
% grid on
set(gca,'FontSize',14)
StandardFigurePBoC([],gca)
box on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
set(gca,'YTick',-1:.25:1)
flux_rate_fig.InvertHardcopy = 'off';
saveas(flux_rate_fig,[FigPath 'rate_constraints.png'])
saveas(flux_rate_fig,[FigPath 'rate_constraints.pdf'])

%% %%%%%%%%%% Plot fractional occupancy for states 3 and 2 vs sharpness %%%
sharpness_norm_index = find(strcmp(metric_names,'SharpnessNormed'));

frac_indices = [sharpness_norm_index frac_index];
metric_one_name = metric_names{frac_indices(1)} ;
metric_two_name = metric_names{frac_indices(2)} ;

% generate load names
load_name_eq = ['param_sweep_results_' metric_one_name '_' metric_two_name ...
             '_eq1.mat'];
         
load_name_neq = ['param_sweep_results_' metric_one_name '_' metric_two_name ...
              '_eq0.mat']; 

% load data
load([DataPath load_name_eq]);
load([DataPath load_name_neq]);  

% extract parameter vectors
metric_array_eq = sim_struct_frac_eq.metric_array;        
mvec1_frac_v_sharp_eq = reshape(metric_array_eq(:,frac_indices(1),:),1,[]);
mvec2_frac_v_sharp_eq = reshape(metric_array_eq(:,frac_indices(2),:),1,[]);
% norm_vec_eq = reshape(metric_array_eq(:,2,:),1,[]);
% norm_vec_eq = norm_vec_eq.*(1-norm_vec_eq);


metric_array_neq = sim_struct_frac_neq.metric_array;        
mvec1_frac_v_sharp_neq = reshape(metric_array_neq(:,frac_indices(1),:),1,[]);
mvec2_frac_v_sharp_neq = reshape(metric_array_neq(:,frac_indices(2),:),1,[]);
% norm_vec_neq = reshape(metric_array_eq(:,2,:),1,[]);
% norm_vec_neq = norm_vec_eq.*(1-norm_vec_eq);

plot_indices_eq = randsample(find(mvec1_frac_v_sharp_eq>=0),2*n_plot,false);
plot_indices_neq = randsample(find(mvec1_frac_v_sharp_neq>=0),2*n_plot,false);            
            
            
frac_figure = figure;

cmap = brewermap(9,'Set2');

hold on
scatter(mvec2_frac_v_sharp_neq(plot_indices_neq),mvec1_frac_v_sharp_neq(plot_indices_neq),...
  markerSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(2,:)); 
scatter(mvec2_frac_v_sharp_eq(plot_indices_eq),mvec1_frac_v_sharp_eq(plot_indices_eq),...
  markerSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(3,:)); 

StandardFigurePBoC([],gca)

xlabel('p_1 + p_3')
ylabel('sharpness (normalized)')
set(gca,'FontSize',14)

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
set(gca,'XTick',0:.2:1)
set(gca,'YTick',0:.25:2)
frac_figure.InvertHardcopy = 'off';
%%
saveas(frac_figure,[FigPath 'sharpness_bounds_2state_frac.png'])
saveas(frac_figure,[FigPath 'sharpness_bounds_2state_frac.pdf'])

% figure;
% scatter(piVecEq(:,1),piVecEq(:,3),20,sharpness_eq_norm)




%% set sim parameters
close all
metric_one_name = 'Flux';
metric_two_name = 'Sharpness';
metric_indices_flux = [find(strcmp(metric_names,metric_one_name)) find(strcmp(metric_names,metric_two_name))];

% generate load names
load_name_flux = ['param_sweep_results_' metric_one_name '_' metric_two_name '_eq0.mat'];

load([DataPath load_name_flux]);
metric_array_flux = sim_struct_flux.metric_array;        
mvec1_flux = reshape(metric_array_flux(:,metric_indices_flux(1),:),1,[]);
mvec2_flux = reshape(metric_array_flux(:,metric_indices_flux(2),:),1,[]);

flux_bins = linspace(0,min(mvec1_flux));
% generate prediction
sharp_pd = 0.25*(1 + (1-exp(flux_bins))./(1+exp(flux_bins))); %NL: this is ad hoc...not sure where this factor of 10 comes from

% make figure
flux_sharp_fig = figure;
cmap2 = brewermap(9,'set2');
hold on
scatter(mvec1_flux,mvec2_flux,markerSize,'MarkerfaceColor',cmap2(5,:),'MarkerEdgeAlpha',0,'MarkerFaceAlpha',markerAlpha);
p1 = plot(flux_bins,sharp_pd,'-','Color','black','LineWidth',2);
p = plot(0,0);

% xlim([0 1.2]);
ylim([0.2 0.55])
xlabel('energy dissipation per cycle (ATP units)')
ylabel('sensitivity (dp_{on}/dc)')
% grid on
set(gca,'FontSize',14)
StandardFigurePBoC(p,gca)
box on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

flux_sharp_fig.InvertHardcopy = 'off';

% saveas(flux_sharp_fig,[FigPath 'sharpness_vs_flux_theory_plot.png'])
% saveas(flux_sharp_fig,[FigPath 'sharpness_vs_flux_theory_plot.pdf'])

%% 

% opt_indices_neq = find(shNeq >= pdNeq.*(1-pdNeq));

%% equilibrium binding plots
ec_vec = linspace(0,10);

piArray = NaN(length(ec_vec),length(ec_vec));
sArray = NaN(length(ec_vec),length(ec_vec));

for c = 1:length(ec_vec)
  for m = 1:length(ec_vec) 
    [piArray(c,m),sArray(c,m)] = eq_sharpness_calculations(ec_vec(c),ec_vec(m));
  end
end