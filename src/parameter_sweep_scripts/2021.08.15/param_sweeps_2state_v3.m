% script to perform simple 2 state parameter sweeps
clear 
close all
addpath(genpath('../utilities/'))

% define save path
OutPath = '../../out/bivariate_parameter_sweeps/';
mkdir(OutPath);
FigPath = '../../fig/two_state_plots/';
mkdir(FigPath)

% Define simulation hypparameters
n_points_unique = 50;
n_points = n_points_unique^2;

% define system parameters
c0 = 1;
c1 = 1.1;
c_true = 1.05;
err_tolerance = .1;
K = log((1-err_tolerance)/err_tolerance);
err_vec = logspace(-3,log10(.25),n_points);
K_vec = log((1-err_vec)./err_vec);
cycle_time = 1e4;

rate_bounds = [-4 4];
koff_u = logspace(rate_bounds(1),rate_bounds(2),n_points_unique); % unbinding range
kon_u = logspace(rate_bounds(1),rate_bounds(2),n_points_unique);%kon_vec = koff_vec./(cycle_time*koff_vec-1)/c_true;%logspace(rate_bounds(1),1,n_points_unique);

koff_vec = ones(1,n_points_unique^2);%repelem(koff_u,n_points_unique);
kon_vec = ones(1,n_points_unique^2);%kon_vec = repmat(kon_u,1,n_points_unique);
cycle_time_vec = (c_true*kon_vec+koff_vec)./(c_true*kon_vec.*koff_vec);

[affinity_vec, sort_order] = sort(kon_vec./koff_vec);
% initialize vectors to store output
dwell_time_vec = NaN(1,n_points);
mRNA_drift_vec = NaN(1,n_points);
micro_drift_vec = NaN(1,n_points);
mRNA_tau_vec = NaN(1,n_points);
micro_tau_vec = NaN(1,n_points);

% conduct sweep
for n = 1:n_points
  kon = kon_vec(n);
  koff = koff_vec(n);
  
  % calculate drift
  [~, ~, ~, ~, ~, ~, V_mRNA_1, ~, V_micro_1, ~, tau_micro_1, tau_mRNA_1]...
                          = calculateDriftStats2State(kon,koff,c0,c1,c1,K_vec(n));
                        
  [~, ~, ~, ~, ~, ~, V_mRNA_0, ~, V_micro_0, ~, tau_micro_0, tau_mRNA_0]...
                          = calculateDriftStats2State(kon,koff,c0,c1,c0,K_vec(n));                        
  
  % add to vectors
  mRNA_drift_vec(n) = .5*V_mRNA_1 + abs(.5*V_mRNA_0);
  micro_drift_vec(n) = .5*V_micro_1 + abs(.5*V_micro_0);
  mRNA_tau_vec(n) = .5*tau_mRNA_1 + .5*tau_mRNA_0;
  micro_tau_vec(n) = .5*tau_micro_1 + .5*tau_micro_0;
end


%% plot cycle time-normalized results

rate_ratio_micro_fig = figure;
cmap = brewermap(9,'Paired');
hold on

% scatter(affinity_vec, micro_drift_vec,'MarkerFaceColor',cmap(4,:),'MarkerEdgeColor','k','MarkerEdgeAlpha',0,'MarkerFaceAlpha',1);%,'LineWidth',2);
% scatter(affinity_vec, mRNA_drift_vec,'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k','MarkerEdgeAlpha',0,'MarkerFaceAlpha',1);%'LineWidth',2);

% plot(affinity_vec, micro_drift_vec(sort_order).*cycle_time_vec(sort_order),'Color',cmap(4,:),'LineWidth',3);
% plot(affinity_vec, mRNA_drift_vec(sort_order).*cycle_time_vec(sort_order),'Color',cmap(3,:),'LineWidth',3);

bar([1],[micro_drift_vec(1).*cycle_time_vec(sort_order(1))],'FaceColor',cmap(4,:))
bar([2],[mRNA_drift_vec(1).*cycle_time_vec(sort_order(1))],'FaceColor',cmap(3,:))
legend('microscopic switching','mRNA only','Location','northeast')

% set(xl ,'Interpreter','latex','FontName','Arial')
ylabel('information per cycle')
grid on
set(gca,'FontSize',14)
set(gca,'xtick',1:2,'xticklabel',{'binding','mRNA'})

set(gca,...
        'Color',[228,221,209]/255) 
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
% ax.XAxis(1).FontSize = 11;
xlim([.5 2.5])

rate_ratio_micro_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(rate_ratio_micro_fig,[FigPath 'two_state_dL_norm.png'])
saveas(rate_ratio_micro_fig,[FigPath 'two_state_dL_norm.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  plot normalized decision time

decision_time_fig = figure;
hold on

% scatter(affinity_vec, micro_tau_vec(si)./cycle_time_vec(si),'MarkerFaceColor',cmap(6,:),'MarkerEdgeColor','k','MarkerEdgeAlpha',0,'MarkerFaceAlpha',1);%,'LineWidth',2);
% scatter(affinity_vec, mRNA_tau_vec(si)./cycle_time_vec(si),'MarkerFaceColor',cmap(5,:),'MarkerEdgeColor','k','MarkerEdgeAlpha',0,'MarkerFaceAlpha',1);%'LineWidth',2);
plot(err_vec, micro_tau_vec(sort_order)./cycle_time_vec(sort_order),'Color',cmap(6,:),'LineWidth',3);
plot(err_vec, mRNA_tau_vec(sort_order)./cycle_time_vec(sort_order),'Color',cmap(5,:),'LineWidth',3);

legend('microscopic switching','mRNA only','Location','northwest')

xl = xlabel('error rate (\epsilon)');
% set(xl ,'Interpreter','latex','FontName','Arial')
ylabel('decision time (number of cycles)')
grid on
set(gca,'FontSize',14)
set(gca,'Xscale','log')
set(gca, 'xdir', 'reverse' )

set(gca,'Color',[228,221,209]/255) 

ax = gca;
% ax.XAxis.MinorTickValues = 10.^(-5:.1:5);
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

% xlim([10^-5 10^1])
ax.XAxis(1).FontSize = 11;
decision_time_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(decision_time_fig,[FigPath 'two_state_decision_time_norm.png'])
saveas(decision_time_fig,[FigPath 'two_state_decision_time_norm.pdf'])