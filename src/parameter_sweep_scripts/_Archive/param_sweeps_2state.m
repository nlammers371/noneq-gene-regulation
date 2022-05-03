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
n_points = 1e3;

% define system parameters
c0 = 1;
c1 = 1.1;

rate_bounds = [-4 4];
tau_cycle = 100;
koff_vec = 1./linspace(0,tau_cycle,n_points); % unbinding range
kon_vec = 1./(tau_cycle-1./koff_vec);

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
                          = calculateDriftStats2State(kon,koff,c0,c1,c1,1);
                        
  [~, ~, ~, ~, ~, ~, V_mRNA_0, ~, V_micro_0, ~, tau_micro_0, tau_mRNA_0]...
                          = calculateDriftStats2State(kon,koff,c0,c1,c0,1);                        
  
  % add to vectors
  mRNA_drift_vec(n) = .5*V_mRNA_1 + abs(.5*V_mRNA_0);
  micro_drift_vec(n) = .5*V_micro_1 + abs(.5*V_micro_0);
  mRNA_tau_vec(n) = .5*tau_mRNA_1 + .5*tau_mRNA_0;
  micro_tau_vec(n) = .5*tau_micro_1 + .5*tau_micro_0;
end

%%
occupancy_vec = .5*c0*kon_vec./(c0*kon_vec+koff_vec) + .5*c1*kon_vec./(c1*kon_vec+koff_vec);
close all

rate_ratio_fig = figure;
cmap = brewermap(9,'Paired');
hold on

scatter(occupancy_vec,micro_drift_vec,'MarkerFaceColor',cmap(4,:),'MarkerEdgeColor','k','MarkerEdgeAlpha',.1);%,'LineWidth',2);
scatter(occupancy_vec,mRNA_drift_vec,'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k','MarkerEdgeAlpha',.05);%'LineWidth',2);

set(gca,'Yscale','log')
legend('microscopic switching','mRNA only','Location','southeast')

StandardFigurePBoC([],gca)

xlabel('fraction of time bound')
ylabel('information rate (dL/ dt)')
set(gca,'FontSize',14)

set(gca,'YTick',fliplr([1e-2 1e-3 1e-4 1e-5 1e-6]))
set(gca,'XTick',0:.25:1)
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

rate_ratio_fig.InvertHardcopy = 'off';

saveas(rate_ratio_fig,[FigPath 'two_state_pon_vs_dL.png'])
saveas(rate_ratio_fig,[FigPath 'two_state_pon_vs_dL.pdf'])


decision_time_fig = figure;
hold on

scatter(occupancy_vec,micro_tau_vec/60,'MarkerFaceColor',cmap(6,:),'MarkerEdgeColor','k','MarkerEdgeAlpha',.1);%,'LineWidth',2);
scatter(occupancy_vec,mRNA_tau_vec/60,'MarkerFaceColor',cmap(5,:),'MarkerEdgeColor','k','MarkerEdgeAlpha',.05);%'LineWidth',2);

set(gca,'Yscale','log')
legend('microscopic switching','mRNA only','Location','northeast')

StandardFigurePBoC([],gca)

xlabel('fraction of time bound')
ylabel('decision time (minutes)')
set(gca,'FontSize',14)

set(gca,'YTick',[1e0 1e1 1e2 1e3 1e4])
set(gca,'XTick',0:.25:1)
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

decision_time_fig.InvertHardcopy = 'off';

saveas(decision_time_fig,[FigPath 'two_state_pon_vs_decision_time.png'])
saveas(decision_time_fig,[FigPath 'two_state_pon_vs_decision_time.pdf'])