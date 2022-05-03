% Plot results of sharpness parameter sweeps
clear 
close all
addpath(genpath('../utilities/'))

% %%%%%%%%%%%%%%%%  set relative read and write paths %%%%%%%%%%%%%%%%%%%%
DataPath = '../../out/bivariate_parameter_sweeps/';
FigPath = '../../fig/info_plots/';
mkdir(FigPath);

% %%%%%%%%%%%%%%%%  get metric names %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,metric_names] = calculateMetrics_v4([]);
flux_index = find(strcmp(metric_names,'Flux'));
rate_index = find(strcmp(metric_names,'Production Rate'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
precision_index = find(strcmp(metric_names,'Precision'));
decision_rate_index = find(strcmp(metric_names,'DecisionRateNorm'));
decision_time_index = find(strcmp(metric_names,'DecisionTimeNorm'));
affinity_index = find(strcmp(metric_names,'Affinity'));
% %%%%%%%%%%%%%%%%  Set plot parameters and constants %%%%%%%%%%%%%%%%%%%%

atp = 50 * 1000 / 6.022e23; % in joules
kbT = 300*1.38e-23;


n_plot = 3e3; % number of points to plot
markerAlpha = 0.5; % marker transparency
markerSize = 25; % marker size

%% %%%%%%%%%%%%%%%%  Plot Affinity vs info rate %%%%%%%%%%%%%%%%%%%%
rng(124)

% set sim parameters
info_indices = [affinity_index decision_rate_index];
metric_one_name = metric_names{info_indices(1)} ;
metric_two_name = metric_names{info_indices(2)} ;

% generate load names
load_name_eq = ['param_sweep_results_' metric_one_name '_' metric_two_name ...
             '_eq1.mat'];
         
load_name_neq = ['param_sweep_results_' metric_one_name '_' metric_two_name ...
              '_eq0.mat'];    
         
% load data
load([DataPath load_name_eq]);
load([DataPath load_name_neq]);    


% extract parameter vectors
metric_array_eq = vertcat(sim_struct_eq.metric_array);        
mvec1_rate_v_aff_eq = metric_array_eq(:,info_indices(1),:);
mvec2_rate_v_aff_eq = metric_array_eq(:,info_indices(2),:);

metric_array_neq = vertcat(sim_struct_neq.metric_array);        
mvec1_rate_v_aff_neq = metric_array_neq(:,info_indices(1),:);
mvec2_rate_v_aff_neq = metric_array_neq(:,info_indices(2),:);

plot_indices_eq = randsample(find(~isnan(mvec1_rate_v_aff_eq)),n_plot,false);
nan_ft_eq = find(~isnan(mvec1_rate_v_aff_eq));
plot_indices_neq = randsample(find(~isnan(mvec1_rate_v_aff_neq)),n_plot,false);
nan_ft_neq = find(~isnan(mvec2_rate_v_aff_eq));
close all

aff_vs_rate = figure;
hold on    
cmap = brewermap(9,'Set2');

neq_boundary_rate = boundary(mvec1_rate_v_aff_neq(nan_ft_neq), mvec2_rate_v_aff_neq(nan_ft_neq),.9);
eq_boundary_rate = boundary(mvec1_rate_v_aff_eq(nan_ft_eq), mvec2_rate_v_aff_eq(nan_ft_eq),.75);


s2 = scatter(0,0,'MarkerFaceColor',cmap(2,:),'MarkerEdgeAlpha',1,'MarkerEdgeColor','k');
scatter(exp(mvec1_rate_v_aff_neq(plot_indices_neq)),mvec2_rate_v_aff_neq(plot_indices_neq),...
  markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(2,:));                

s1 = scatter(0,0,'MarkerFaceColor',cmap(3,:),'MarkerEdgeAlpha',1,'MarkerEdgeColor','k');
scatter(exp(mvec1_rate_v_aff_eq(plot_indices_eq)),mvec2_rate_v_aff_eq(plot_indices_eq),...
  markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(3,:));  

p2 = patch(exp(mvec1_rate_v_aff_neq(nan_ft_neq(neq_boundary_rate))), mvec2_rate_v_aff_neq(nan_ft_neq(neq_boundary_rate)),...
  cmap(2,:),'EdgeColor','k','LineWidth',1.5);
p2.FaceAlpha = 0.1;

p1 = patch(exp(mvec1_rate_v_aff_eq(nan_ft_eq(eq_boundary_rate))), mvec2_rate_v_aff_eq(nan_ft_eq(eq_boundary_rate)),...
  cmap(3,:),'EdgeColor','k','LineWidth',1.5);
p1.FaceAlpha = 0.1;
                

box on    

xl = xlabel('TF affinity ($\frac{k_{on}}{k_{off}}$)');
set(xl ,'Interpreter','latex')
yl = ylabel('information per cycle');
grid on
set(gca,'FontSize',14)
set(gca,'Xscale','log')
set(gca,'Color',[228,221,209]/255) 
ax = gca;
ax.YAxis(1).Color = 'k';

% xlim([1e-6 1e6])
ylim([0 0.01])

aff_vs_rate.InvertHardcopy = 'off';
set(gcf,'color','w');

legend([s1 s2],'equilibrium','non-equilibrium','Color','w','Location','southeast')
saveas(aff_vs_rate,[FigPath metric_one_name '_' metric_two_name '_neq.png'])
saveas(aff_vs_rate,[FigPath metric_one_name '_' metric_two_name '_neq.pdf'])


%% %%%%%%%%%%%%%%%%  Plot Affinity vs decision time %%%%%%%%%%%%%%%%%%%%

% set sim parameters
info_indices = [affinity_index decision_time_index];
metric_one_name = metric_names{info_indices(1)} ;
metric_two_name = metric_names{info_indices(2)} ;

% generate load names
load_name_eq = ['param_sweep_results_' metric_one_name '_' metric_two_name ...
             '_eq1.mat'];
         
load_name_neq = ['param_sweep_results_' metric_one_name '_' metric_two_name ...
              '_eq0.mat'];    
         
% load data
load([DataPath load_name_eq]);
load([DataPath load_name_neq]);    

% extract parameter vectors
metric_array_eq = vertcat(sim_struct_eq.metric_array);        
mvec1_rate_v_aff_eq = metric_array_eq(:,info_indices(1),:);
mvec2_rate_v_aff_eq = 1./metric_array_eq(:,info_indices(2),:);

metric_array_neq = vertcat(sim_struct_neq.metric_array);        
mvec1_rate_v_aff_neq = metric_array_neq(:,info_indices(1),:);
mvec2_rate_v_aff_neq = 1./metric_array_neq(:,info_indices(2),:);

plot_indices_eq = randsample(find(~isnan(mvec1_rate_v_aff_eq)),n_plot,false);
nan_ft_eq = find(~isnan(mvec1_rate_v_aff_eq)&~isnan(mvec2_rate_v_aff_eq));
plot_indices_neq = randsample(find(~isnan(mvec1_rate_v_aff_neq)),n_plot,false);
nan_ft_neq = find(~isnan(mvec1_rate_v_aff_neq)&~isnan(mvec2_rate_v_aff_neq));
close all

aff_vs_rate = figure;
hold on    
cmap = brewermap(9,'Set2');

neq_boundary_rate = boundary(mvec1_rate_v_aff_neq(nan_ft_neq), mvec2_rate_v_aff_neq(nan_ft_neq),.85);
eq_boundary_rate = boundary(mvec1_rate_v_aff_eq(nan_ft_eq), mvec2_rate_v_aff_eq(nan_ft_eq),.75);


s2 = scatter(0,0,'MarkerFaceColor',cmap(2,:),'MarkerEdgeAlpha',1,'MarkerEdgeColor','k');
scatter(exp(mvec1_rate_v_aff_neq(plot_indices_neq)),mvec2_rate_v_aff_neq(plot_indices_neq),...
  markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(2,:));                

s1 = scatter(0,0,'MarkerFaceColor',cmap(3,:),'MarkerEdgeAlpha',1,'MarkerEdgeColor','k');
scatter(exp(mvec1_rate_v_aff_eq(plot_indices_eq)),mvec2_rate_v_aff_eq(plot_indices_eq),...
  markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(3,:));  

p2 = patch(exp(mvec1_rate_v_aff_neq(nan_ft_neq(neq_boundary_rate))), mvec2_rate_v_aff_neq(nan_ft_neq(neq_boundary_rate)),...
  cmap(2,:),'EdgeColor','k','LineWidth',1.5);
p2.FaceAlpha = 0.1;

p1 = patch(exp(mvec1_rate_v_aff_eq(nan_ft_eq(eq_boundary_rate))), mvec2_rate_v_aff_eq(nan_ft_eq(eq_boundary_rate)),...
  cmap(3,:),'EdgeColor','k','LineWidth',1.5);
p1.FaceAlpha = 0.1;
                

box on    

xl = xlabel('TF affinity ($\frac{k_{on}}{k_{off}}$)');
set(xl ,'Interpreter','latex')
yl = ylabel('decision time (number of cycles)');
grid on
set(gca,'FontSize',14)
set(gca,'Xscale','log')
set(gca,'Color',[228,221,209]/255) 
ax = gca;
ax.YAxis(1).Color = 'k';

% xlim([1e-6 1e6])
ylim([0 2e3])

aff_vs_rate.InvertHardcopy = 'off';
set(gcf,'color','w');

legend([s1 s2],'equilibrium','non-equilibrium','Color','w','Location','northeast')
saveas(aff_vs_rate,[FigPath metric_one_name '_' metric_two_name '_neq.png'])
saveas(aff_vs_rate,[FigPath metric_one_name '_' metric_two_name '_neq.pdf'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%% sharpness vs noise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
info_indices2 = [sharpness_index precision_index];

metric_one_name = metric_names{info_indices2(1)} ;
metric_two_name = metric_names{info_indices2(2)} ;


% generate load names
load_name_eq = ['param_sweep_results_' metric_one_name '_' metric_two_name ...
             '_eq1.mat'];
         
load_name_neq = ['param_sweep_results_' metric_one_name '_' metric_two_name ...
              '_eq0.mat'];    
         
% load data
load([DataPath load_name_eq]);
load([DataPath load_name_neq]);

metric_array_eq = vertcat(sim_struct_eq.metric_array);      
mvec1_sharp_v_noise_eq = metric_array_eq(:,info_indices2(1),:);
mvec2_sharp_v_noise_eq = 1./metric_array_eq(:,info_indices2(2),:);

metric_array_neq = vertcat(sim_struct_neq.metric_array);      
mvec1_sharp_v_noise_neq = metric_array_neq(:,info_indices2(1),:);
mvec2_sharp_v_noise_neq = 1./metric_array_neq(:,info_indices2(2),:);

plot_indices_eq = randsample(find(mvec1_sharp_v_noise_eq>=0),n_plot,false);
plot_indices_neq = randsample(find(mvec1_sharp_v_noise_neq>=0),n_plot,false);
close all


sharpness_vs_precision = figure;
hold on    
cmap = brewermap(9,'Set2');


s1 = scatter(mvec1_sharp_v_noise_neq(plot_indices_neq),mvec2_sharp_v_noise_neq(plot_indices_neq),...
  markerSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(2,:));                

s2 = scatter(mvec1_sharp_v_noise_eq(plot_indices_eq),mvec2_sharp_v_noise_eq(plot_indices_eq),...
  markerSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(3,:));                

box on    

xl = xlabel('sharpness (dr / dc)');
yl = ylabel('intrinsic noise (\sigma_{r})');
set(yl ,'Interpreter','latex')
grid on
set(gca,'FontSize',14)
%   set(gca,'Xscale','log')
set(gca,'Color',[228,221,209]/255) 
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
ax.XAxis(1).FontSize = 11;

xlim([0 0.5])
ylim([0 10])
set(gca,'Yscale','log')

sharpness_vs_precision.InvertHardcopy = 'off';
set(gcf,'color','w');

legend([s2 s1],'equilibrium','non-equilibrium','Color','w','Location','northeast')
saveas(sharpness_vs_precision,[FigPath metric_one_name '_' metric_two_name '_neq.png'])
saveas(sharpness_vs_precision,[FigPath metric_one_name '_' metric_two_name '_neq.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%% sharpness vs noise (HM) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
info_indices2 = [sharpness_index precision_index];

metric_one_name = metric_names{info_indices2(1)} ;
metric_two_name = metric_names{info_indices2(2)} ;


% generate load names
load_name_eq = ['param_sweep_results_' metric_one_name '_' metric_two_name ...
             '_half_eq1.mat'];
         
load_name_neq = ['param_sweep_results_' metric_one_name '_' metric_two_name ...
              '_half_eq0.mat'];    
         
% load data
load([DataPath load_name_eq]);
load([DataPath load_name_neq]);

metric_array_eq = vertcat(sim_struct_eq.metric_array);      
mvec1_sharp_v_noise_eq = metric_array_eq(:,info_indices2(1),:);
mvec2_sharp_v_noise_eq = metric_array_eq(:,info_indices2(2),:);
pd_vec_eq = metric_array_eq(:,rate_index,:);

metric_array_neq = vertcat(sim_struct_neq.metric_array);      
mvec1_sharp_v_noise_neq = metric_array_neq(:,info_indices2(1),:);
mvec2_sharp_v_noise_neq = metric_array_neq(:,info_indices2(2),:);
pd_vec_neq = metric_array_neq(:,rate_index,:);

plot_indices_eq = randsample(find(mvec1_sharp_v_noise_eq>=0&round(pd_vec_eq,1)==.5),n_plot,false);
plot_indices_neq = randsample(find(mvec1_sharp_v_noise_neq>=0&round(pd_vec_neq,1)==.5),n_plot,false);
close all


sharpness_vs_precision = figure;
hold on    
cmap = brewermap(9,'Set2');


s1 = scatter(mvec1_sharp_v_noise_neq(plot_indices_neq),mvec2_sharp_v_noise_neq(plot_indices_neq),...
  markerSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(2,:));                

s2 = scatter(mvec1_sharp_v_noise_eq(plot_indices_eq),mvec2_sharp_v_noise_eq(plot_indices_eq),...
  markerSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(3,:));                

box on    

xl = xlabel('sharpness (dr / dc)');
yl = ylabel('intrinsic noise (\sigma_{r})');
set(yl ,'Interpreter','latex')
grid on
set(gca,'FontSize',14)
%   set(gca,'Xscale','log')
set(gca,'Color',[228,221,209]/255) 
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
ax.XAxis(1).FontSize = 11;

xlim([0 0.5])
ylim([0 10])
set(gca,'Yscale','log')

sharpness_vs_precision.InvertHardcopy = 'off';
set(gcf,'color','w');

legend([s2 s1],'equilibrium','non-equilibrium','Color','w','Location','northeast')
saveas(sharpness_vs_precision,[FigPath metric_one_name '_' metric_two_name '_HM.png'])
saveas(sharpness_vs_precision,[FigPath metric_one_name '_' metric_two_name '_HM.pdf'])
  
%% %%%%%%%%%%%%%%%%%%%% Information vs flux %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load_name_flux = ['param_sweep_results_' metric_names{flux_index} '_' ...
    metric_names{decision_rate_index} '_eq0.mat'];

        
% load data
load([DataPath load_name_flux]);

% extract vectors
metric_array_neq = vertcat(sim_struct_neq.metric_array);      
flux_vec = metric_array_neq(:,flux_index,:);
info_vec = metric_array_neq(:,decision_rate_index,:);
% nan_ft = find(flux_vec>=0);
pd_vec_neq = metric_array_neq(:,rate_index,:);

% randomly draw indices to show
plot_indices_neq = randsample(find(flux_vec>=0),n_plot,false);

% find boundary
neq_boundary = boundary(flux_vec(plot_indices_neq), info_vec(plot_indices_neq),.5);

close all

energy_v_info = figure;
hold on    
cmap = brewermap(9,'Set2');

scatter(flux_vec(plot_indices_neq), info_vec(plot_indices_neq),...
          markerSize,'MarkerEdgeAlpha',.2,'MarkerEdgeColor','k','MarkerFaceAlpha',...
          markerAlpha, 'MarkerFaceColor',cmap(5,:));                

p1 = patch(flux_vec(plot_indices_neq(neq_boundary)), info_vec(plot_indices_neq(neq_boundary)),...
          cmap(5,:),'EdgeColor','k','LineWidth',1.5);      
p1.FaceAlpha = 0.1;

box on    

xl = xlabel('dissipation per cycle (k_bT)');
yl = ylabel('information per cycle');

grid on
set(gca,'FontSize',14)
%   set(gca,'Xscale','log')
set(gca,'Color',[228,221,209]/255) 
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
ax.XAxis(1).FontSize = 11;

xlim([0 30])
% ylim([0 10])

energy_v_info.InvertHardcopy = 'off';
set(gcf,'color','w');


saveas(energy_v_info,[FigPath 'Flux_DecisionRateNorn.png'])
saveas(energy_v_info,[FigPath 'Flux_DecisionRateNorn.pdf'])  
%  
% %% %%%%%%%%%% Plot on and off rate c dependencies %%%%%%%%%%%%%%%%%%%%%%%%%
% close all
% 
% % get variable indices
% kon_index = find(strcmp(metric_names,'KonEff'));
% koff_index = find(strcmp(metric_names,'KoffEff'));
% dkon_index = find(strcmp(metric_names,'dKondC'));
% dkoff_index = find(strcmp(metric_names,'dKoffdC'));
% dkdc_indices = [dkon_index dkoff_index];
% % name of data sets
% save_name_dkdc = ['param_sweep_results_' metric_names{dkdc_indices(1)} '_' ...
%     metric_names{dkdc_indices(2)}];
% 
% % load
% eq_flags = [0 1];
% dkdc_struct = struct;
% for e = 1:length(eq_flags)
%   load([DataPath save_name_dkdc '_eq' num2str(eq_flags(e)) '.mat'],'sim_struct')
%   dkdc_struct(e).sim_struct = sim_struct;
% end
% 
% % off and on rate c dependencies
% koffNeq = reshape(dkdc_struct(1).sim_struct.metric_array(:,koff_index,:),1,[]);
% konNeq = reshape(dkdc_struct(1).sim_struct.metric_array(:,kon_index,:),1,[]);
% dKoffdCNeq = reshape(dkdc_struct(1).sim_struct.metric_array(:,dkoff_index,:),1,[]);%./koffNeq;
% dKondCNeq = reshape(dkdc_struct(1).sim_struct.metric_array(:,dkon_index,:),1,[]);%./konNeq;
% 
% pdRateEq = reshape(dkdc_struct(2).sim_struct.metric_array(:,2,:),1,[]);
% sharpnessEq = reshape(dkdc_struct(2).sim_struct.metric_array(:,1,:),1,[])./(pdRateEq.*(1-pdRateEq));
% fracEq = reshape(dkdc_struct(2).sim_struct.metric_array(:,frac_index,:),1,[]);
% 
% koffEq = reshape(dkdc_struct(2).sim_struct.metric_array(:,koff_index,:),1,[]);
% konEq = reshape(dkdc_struct(2).sim_struct.metric_array(:,kon_index,:),1,[]);
% dKoffdCEq = reshape(dkdc_struct(2).sim_struct.metric_array(:,dkoff_index,:),1,[]);%./koffEq;
% dKondCEq = reshape(dkdc_struct(2).sim_struct.metric_array(:,dkon_index,:),1,[]);%./konEq;
% 
% plot_indices_eq = randsample(find(sharpnessEq>=0.0),min([n_plot length(find(sharpnessEq>=0.0))]),false);
% 
% fluxNeq = abs(reshape(dkdc_struct(1).sim_struct.metric_array(:,3,:),1,[]));
% pdRateNeq = reshape(dkdc_struct(1).sim_struct.metric_array(:,2,:),1,[]);
% sharpnessNeq = reshape(dkdc_struct(1).sim_struct.metric_array(:,1,:),1,[])./(pdRateNeq.*(1-pdRateNeq));
% 
% plot_indices_neq = randsample(find(sharpnessNeq>=0.0),min([n_plot length(find(sharpnessNeq>=0.0))]),false);
% 
% 
% % make figure
% flux_rate_fig = figure;
% cmap3 = flipud(brewermap([],'Spectral'));
% hold on
% colormap(cmap3);
% 
% scatter(dKondCNeq(plot_indices_neq),dKoffdCNeq(plot_indices_neq),markerSize,sharpnessNeq(plot_indices_neq),'filled',...
%         'MarkerEdgeAlpha',.1,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha);
% 
% scatter(dKondCEq(plot_indices_eq),dKoffdCEq(plot_indices_eq),markerSize,sharpnessEq(plot_indices_eq),'filled',...
%         'MarkerEdgeAlpha',.1,'MarkerFaceAlpha',markerAlpha,'MarkerEdgeColor','k');
% 
% 
% h = colorbar;
% ylabel(h,'normalized sharpness')
% caxis([0 2])
% xlim([0 2.1]);
% ylim([-1 0]);
% % ylim([-1 0])
% xlabel('d k_{on} / dc')
% ylabel('d k_{off} / dc')
% % grid on
% set(gca,'FontSize',14)
% StandardFigurePBoC([],gca)
% box on
% ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.XAxis(1).Color = 'k';
% set(gca,'YTick',-1:.25:1)
% flux_rate_fig.InvertHardcopy = 'off';
% saveas(flux_rate_fig,[FigPath 'rate_constraints.png'])
% saveas(flux_rate_fig,[FigPath 'rate_constraints.pdf'])
% 
% %% %%%%%%%%%% Plot fractional occupancy for states 3 and 2 vs sharpness %%%
% sharpness_norm_index = find(strcmp(metric_names,'SharpnessNormed'));
% 
% frac_indices = [sharpness_norm_index frac_index];
% metric_one_name = metric_names{frac_indices(1)} ;
% metric_two_name = metric_names{frac_indices(2)} ;
% 
% % generate load names
% load_name_eq = ['param_sweep_results_' metric_one_name '_' metric_two_name ...
%              '_eq1.mat'];
%          
% load_name_neq = ['param_sweep_results_' metric_one_name '_' metric_two_name ...
%               '_eq0.mat']; 
% 
% % load data
% load([DataPath load_name_eq]);
% load([DataPath load_name_neq]);  
% 
% % extract parameter vectors
% metric_array_eq = sim_struct_frac_eq.metric_array;        
% mvec1_frac_v_sharp_eq = reshape(metric_array_eq(:,frac_indices(1),:),1,[]);
% mvec2_frac_v_sharp_eq = reshape(metric_array_eq(:,frac_indices(2),:),1,[]);
% % norm_vec_eq = reshape(metric_array_eq(:,2,:),1,[]);
% % norm_vec_eq = norm_vec_eq.*(1-norm_vec_eq);
% 
% 
% metric_array_neq = sim_struct_frac_neq.metric_array;        
% mvec1_frac_v_sharp_neq = reshape(metric_array_neq(:,frac_indices(1),:),1,[]);
% mvec2_frac_v_sharp_neq = reshape(metric_array_neq(:,frac_indices(2),:),1,[]);
% % norm_vec_neq = reshape(metric_array_eq(:,2,:),1,[]);
% % norm_vec_neq = norm_vec_eq.*(1-norm_vec_eq);
% 
% plot_indices_eq = randsample(find(mvec1_frac_v_sharp_eq>=0),2*n_plot,false);
% plot_indices_neq = randsample(find(mvec1_frac_v_sharp_neq>=0),2*n_plot,false);            
%             
%             
% frac_figure = figure;
% 
% cmap = brewermap(9,'Set2');
% 
% hold on
% scatter(mvec2_frac_v_sharp_neq(plot_indices_neq),mvec1_frac_v_sharp_neq(plot_indices_neq),...
%   markerSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(2,:)); 
% scatter(mvec2_frac_v_sharp_eq(plot_indices_eq),mvec1_frac_v_sharp_eq(plot_indices_eq),...
%   markerSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(3,:)); 
% 
% StandardFigurePBoC([],gca)
% 
% xlabel('p_1 + p_3')
% ylabel('sharpness (normalized)')
% set(gca,'FontSize',14)
% 
% ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.XAxis(1).Color = 'k';
% set(gca,'XTick',0:.2:1)
% set(gca,'YTick',0:.25:2)
% frac_figure.InvertHardcopy = 'off';
% %%
% saveas(frac_figure,[FigPath 'sharpness_bounds_2state_frac.png'])
% saveas(frac_figure,[FigPath 'sharpness_bounds_2state_frac.pdf'])
% 
% % figure;
% % scatter(piVecEq(:,1),piVecEq(:,3),20,sharpness_eq_norm)
% 
% 
% 
% 
% %% set sim parameters
% close all
% metric_one_name = 'Flux';
% metric_two_name = 'Sharpness';
% metric_indices_flux = [find(strcmp(metric_names,metric_one_name)) find(strcmp(metric_names,metric_two_name))];
% 
% % generate load names
% load_name_flux = ['param_sweep_results_' metric_one_name '_' metric_two_name '_eq0.mat'];
% 
% load([DataPath load_name_flux]);
% metric_array_flux = sim_struct_flux.metric_array;        
% mvec1_flux = reshape(metric_array_flux(:,metric_indices_flux(1),:),1,[]);
% mvec2_flux = reshape(metric_array_flux(:,metric_indices_flux(2),:),1,[]);
% 
% flux_bins = linspace(0,min(mvec1_flux));
% % generate prediction
% sharp_pd = 0.25*(1 + (1-exp(flux_bins))./(1+exp(flux_bins))); %NL: this is ad hoc...not sure where this factor of 10 comes from
% 
% % make figure
% flux_sharp_fig = figure;
% cmap2 = brewermap(9,'set2');
% hold on
% scatter(mvec1_flux,mvec2_flux,markerSize,'MarkerfaceColor',cmap2(5,:),'MarkerEdgeAlpha',0,'MarkerFaceAlpha',markerAlpha);
% p1 = plot(flux_bins,sharp_pd,'-','Color','black','LineWidth',2);
% p = plot(0,0);
% 
% % xlim([0 1.2]);
% ylim([0.2 0.55])
% xlabel('energy dissipation per cycle (ATP units)')
% ylabel('sensitivity (dp_{on}/dc)')
% % grid on
% set(gca,'FontSize',14)
% StandardFigurePBoC(p,gca)
% box on
% ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.XAxis(1).Color = 'k';
% 
% flux_sharp_fig.InvertHardcopy = 'off';
% 
% % saveas(flux_sharp_fig,[FigPath 'sharpness_vs_flux_theory_plot.png'])
% % saveas(flux_sharp_fig,[FigPath 'sharpness_vs_flux_theory_plot.pdf'])
% 
% %% 
% 
% % opt_indices_neq = find(shNeq >= pdNeq.*(1-pdNeq));
% 
% %% equilibrium binding plots
% ec_vec = linspace(0,10);
% 
% piArray = NaN(length(ec_vec),length(ec_vec));
% sArray = NaN(length(ec_vec),length(ec_vec));
% 
% for c = 1:length(ec_vec)
%   for m = 1:length(ec_vec) 
%     [piArray(c,m),sArray(c,m)] = eq_sharpness_calculations(ec_vec(c),ec_vec(m));
%   end
% end