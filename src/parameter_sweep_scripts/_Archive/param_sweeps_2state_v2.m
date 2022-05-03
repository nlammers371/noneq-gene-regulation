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
c_true = 1;

rate_bounds = [-1 4];
cycle_time = 60; % set to 1 minute
koff_u = logspace(rate_bounds(1),rate_bounds(2),n_points_unique^2); % unbinding range
kon_u = koff_u./(cycle_time*koff_u-1);%logspace(rate_bounds(1),1,n_points_unique);

koff_vec = koff_u;%repmat(koff_u,1,n_points_unique);
kon_vec = kon_u;%repelem(kon_u,n_points_unique);

% initialize vectors to store output
dwell_time_vec = NaN(1,n_points);
mRNA_noise_vec = NaN(1,n_points);
micro_noise_vec = NaN(1,n_points);
mRNA_tau_vec = NaN(1,n_points);
micro_tau_vec = NaN(1,n_points);

% conduct sweep
for n = 1:n_points
  kon = kon_vec(n);
  koff = koff_vec(n);
  
  % calculate variance
  [r, v_mRNA, v_ML] = calculateVar2State(kon,koff,c_true);
  
  % add to vectors
  mRNA_noise_vec(n) = v_mRNA;% + abs(.5*V_mRNA_0);
  micro_noise_vec(n) = v_ML;% + abs(.5*V_micro_0);  
end

%%
cycle_time_vec = (c_true*kon_vec+koff_vec)./(c_true*kon_vec.*koff_vec);% + .5*c1*kon_vec.*koff_vec./(c1*kon_vec+koff_vec)); 
% find Bcd-like stuff
Bcd_filter = koff_vec>0.5 & koff_vec<.7 ;%& kon_vec <= 1;
yBoundsNoise = [.5e-6 1.1e0];
xBounds = [10^-3.1 10^3.1];
[affinity_vec, si] = sort(kon_vec./koff_vec);
close all

% plot cycle time-normalized results

rate_ratio_micro_fig = figure;
cmap = brewermap(9,'Paired');
hold on

scatter(affinity_vec, 100*(micro_noise_vec(si)/60),'MarkerFaceColor',...
  cmap(4,:),'MarkerEdgeColor','k','MarkerEdgeAlpha',0,'MarkerFaceAlpha',1);%,'LineWidth',2);
scatter(affinity_vec, 100*(mRNA_noise_vec(si)/60),'MarkerFaceColor',...
  cmap(3,:),'MarkerEdgeColor','k','MarkerEdgeAlpha',0,'MarkerFaceAlpha',1);%'LineWidth',2);

legend('microscopic switching','mRNA only','Location','southwest')

xl = xlabel('TF affinity (k_{on}/k_{off})');
% set(xl ,'Interpreter','latex','FontName','Arial')
yl = ylabel('decision time (minutes)');
% set(yl ,'Interpreter','latex','FontName','Arial')
grid on
set(gca,'FontSize',14)
set(gca,'Xscale','log')
% set(gca,'Yscale','log')
set(gca,'Color',[228,221,209]/255) 
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).FontSize = 11;
ylim([25 250])
% xlim(xBounds)

rate_ratio_micro_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(rate_ratio_micro_fig,[FigPath 'two_state_noise_norm.png'])
saveas(rate_ratio_micro_fig,[FigPath 'two_state_noise_norm.pdf'])
