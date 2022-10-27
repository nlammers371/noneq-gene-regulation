% This script to make figures examining validity of Gaussian noise
% approximation
clear
close all
addpath(genpath(['..' filesep 'utilities']))
% DropboxFolder = ['C:' filesep 'Users' filesep 'nlamm' filesep 'Dropbox (Personal)' filesep 'Nonequilibrium' filesep 'Nick' filesep];
DropboxFolder = [filesep 'Users' filesep 'nick' filesep 'Dropbox (Personal)' filesep 'Nonequilibrium' filesep 'Nick' filesep];
ReadPath = [DropboxFolder  'SweepOutput'  filesep 'appendices' filesep];
FigPath = [DropboxFolder 'manuscript' filesep 'appendices' filesep 'gaussian_approximation' filesep];
mkdir(FigPath)

% load simulation data
load([ReadPath 'ss_struct.mat'])

%% Calculate percent absolute deviation from ground truth
var_dev_array = NaN(length(ss_struct(1).r_mean_vec),length(ss_struct));
mean_dev_array = NaN(length(ss_struct(1).r_mean_vec),length(ss_struct));

var_abs_dev_array = NaN(length(ss_struct(1).r_mean_vec),length(ss_struct));
mean_abs_dev_array = NaN(length(ss_struct(1).r_mean_vec),length(ss_struct));

var_pd_array = NaN(length(ss_struct(1).r_mean_vec),length(ss_struct));
mean_pd_array = NaN(length(ss_struct(1).r_mean_vec),length(ss_struct));

var_sim_array = NaN(length(ss_struct(1).r_mean_vec),length(ss_struct));
mean_sim_array = NaN(length(ss_struct(1).r_mean_vec),length(ss_struct));

for i = 1:length(ss_struct)
    rt = ss_struct(i).r_predicted;
    vt = ss_struct(i).var_predicted;

    mean_dev_array(:,i) = (rt-ss_struct(i).r_mean_vec)/rt;
    var_dev_array(:,i) = (vt-ss_struct(i).r_var_vec)/vt;

    mean_abs_dev_array(:,i) = abs(rt-ss_struct(i).r_mean_vec)/rt;
    var_abs_dev_array(:,i) = abs(vt-ss_struct(i).r_var_vec)/vt;

    % record values for subsequent use
    mean_sim_array(:,i) = ss_struct(i).r_mean_vec;
    var_sim_array(:,i) = ss_struct(i).r_var_vec;

    mean_pd_array(:,i) = rt;
    var_pd_array(:,i) = vt;

end

mean_abs_dev_array = mean_abs_dev_array - mean_abs_dev_array(end,:);
var_abs_dev_array = var_abs_dev_array - var_abs_dev_array(end,:);
%% now look at mean and variance values for gene circuits, grouped by mean rate
n_groups = 10;
r_bins = linspace(0,1,n_groups+1);
r_mean_vec = [ss_struct.r_mean];
r_group_vec = discretize(r_mean_vec,r_bins);
% p_array_sim = p_array';%vertcat(ss_struct.p_vec);

% use bootstrapping to estimate standard error in convergence across
% different groups
t_grid_long = ss_struct(1).time_grid;
nBoots = 100;
ma_boot_array = NaN(length(t_grid_long),n_groups,nBoots);
% m2_boot_array = NaN(length(t_grid_long),n_groups,nBoots);
m_boot_array = NaN(length(t_grid_long),n_groups,nBoots);
va_boot_array = NaN(length(t_grid_long),n_groups,nBoots);
% v2_boot_array = NaN(length(t_grid_long),n_groups,nBoots);
v_boot_array = NaN(length(t_grid_long),n_groups,nBoots);

for n = 1:n_groups
    option_vec = find(r_group_vec==n);
    for b = 1:nBoots
        boot_indices = randsample(option_vec,length(option_vec),true);

        % estimate percent deviation statistics
        ma_boot_array(:,n,b) = mean(mean_abs_dev_array(:,boot_indices),2);
        m_boot_array(:,n,b) = mean(mean_dev_array(:,boot_indices),2);
        va_boot_array(:,n,b) = mean(var_abs_dev_array(:,boot_indices),2);
        v_boot_array(:,n,b) = mean(var_dev_array(:,boot_indices),2);

    end
end

m_mean_array = mean(m_boot_array,3)*100;
m_ste_array = std(m_boot_array,[],3)*100;

ma_mean_array = mean(ma_boot_array,3)*100;
ma_ste_array = std(ma_boot_array,[],3)*100;

v_mean_array = mean(v_boot_array,3)*100;
v_ste_array = std(v_boot_array,[],3)*100;

va_mean_array = mean(va_boot_array,3)*100;
va_ste_array = std(va_boot_array,[],3)*100;


close all
ma_fig = figure;
hold on
cmap = flipud(brewermap(n_groups,'Spectral'));
colormap(cmap); 
for i = 1:n_groups
    errorbar(t_grid_long,ma_mean_array(:,i),ma_ste_array(:,i),'Color',cmap(i,:),'Capsize',0,'LineWidth',1.5)
end  

grid on
xlabel('time (burst cycles)');
ylabel('percent absolute deviation (production rate)')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
h = colorbar;
ylabel(h,'average production rate (r)')
ax = gca;
xlim([1e-1 5e3])
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

ma_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

ma_fig.Renderer='Painters';
set(gca,'xscale','log')
% set(gca,'yscale','log')
ylim([-1 30])
saveas(ma_fig,[FigPath,'ma_trend_log.png'])
saveas(ma_fig,[FigPath,'ma_trend_log.pdf'])

% variance
va_fig = figure;
hold on
cmap = flipud(brewermap(n_groups,'Spectral'));
colormap(cmap); 
for i = 1:n_groups
    errorbar(t_grid_long,va_mean_array(:,i),va_ste_array(:,i),'Color',cmap(i,:),'Capsize',0,'LineWidth',1.5)
end  

grid on
xlabel('time (burst cycles)');
ylabel('percent absolute deviation (variance)')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
h = colorbar;
ylabel(h,'average production rate (r)')
ax = gca;
xlim([1e-1 5e3])
ylim([-5 80])
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

va_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

va_fig.Renderer='Painters';
set(gca,'xscale','log')
% set(gca,'yscale','log')

saveas(va_fig,[FigPath,'va_trend_log.png'])
saveas(va_fig,[FigPath,'va_trend_log.pdf'])


%% make scatter plots showing predicted vs actual for final time point
plot_ind = find(t_grid_long>=32.5,1);
plot_ind_2 = find(t_grid_long>=100,1);
rng(123);

r_scatter_fig = figure;
hold on
cmap = brewermap([],'set2');
plot_options = 1:length(mean_pd_array(plot_ind,:));
plot_1 = randsample(plot_options,length(plot_options),false);
plot_2 = randsample(plot_options,length(plot_options),false);
for p = 1:10:length(plot_options)
    scatter(mean_pd_array(end,plot_1(p:p+10-1)),mean_sim_array(end,plot_1(p:p+10-1)),75,'d','MarkerFaceColor',cmap(8,:),'MarkerEdgeColor','k','MarkerFaceAlpha',1)
    scatter(mean_pd_array(plot_ind,plot_2(p:p+10-1)),mean_sim_array(plot_ind,plot_2(p:p+10-1)),75,'MarkerFaceColor',cmap(5,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.5)
end
plot(linspace(0,1),linspace(0,1),'-','Color','k','LineWidth',2)

% grid on
xlabel('production rate (analytic prediction)');
ylabel('production rate (simulation)')


% set(gca,'Color',[228,221,209]/255) 
legend('32.5 burst cycles', '100 burst cycles','Location','southeast')

set(gca,'FontSize',14)

set(gca,'xtick',0:0.2:1,'ytick',0:0.2:1)
ax = gca;
xlim([0 1])
ylim([0 1])
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
box on
r_scatter_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(r_scatter_fig,[FigPath,'r_scatter.png'])
saveas(r_scatter_fig,[FigPath,'r_scatter.pdf'])


min_v = min(var_pd_array(plot_ind,:));
max_v = max(var_pd_array(plot_ind,:));
%%
v_scatter_fig = figure;
hold on
cmap = brewermap([],'set2');

% scatter(var_pd_array(plot_ind,:),var_sim_array(plot_ind_2,:),75,'d','MarkerFaceColor',cmap(8,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.5)
% scatter(var_pd_array(plot_ind,:),var_sim_array(plot_ind,:),75,'MarkerFaceColor',cmap(5,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.5)
for p = 1:10:length(plot_options)
    scatter(var_pd_array(end,plot_1(p:p+10-1)),var_sim_array(end,plot_1(p:p+10-1)),75,'d','MarkerFaceColor',cmap(8,:),'MarkerEdgeColor','k','MarkerFaceAlpha',1)
    scatter(var_pd_array(plot_ind,plot_2(p:p+10-1)),var_sim_array(plot_ind,plot_2(p:p+10-1)),75,'MarkerFaceColor',cmap(5,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.5)
end
plot(logspace(-4,1),logspace(-4,1),'-.','Color','k','LineWidth',2)

% grid on
xlabel('variance (analytic prediction)');
ylabel('variance (simulation)')
box on
set(gca,'FontSize',14)
% set(gca,'Color',[228,221,209]/255) 

% set(gca,'xtick',0:0.2:1,'ytick',0:0.2:1)

set(gca,'yscale','log')
set(gca,'xscale','log')

ax = gca;
xlim([1e-4 1e1])
ylim([1e-4 1e1])
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

legend('32.5 burst cycles', '5000 burst cycles','Location','southeast')

v_scatter_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(v_scatter_fig,[FigPath,'v_scatter.png'])
saveas(v_scatter_fig,[FigPath,'v_scatter.pdf'])

var_err32 = mean(var_abs_dev_array(plot_ind,:))
mean_err32 = mean(mean_abs_dev_array(plot_ind,:))

var_err_long = mean(mean(var_abs_dev_array(plot_ind_2:end,:)))
var_err_ste_long = std(mean(var_abs_dev_array(plot_ind_2:end,:)))
mean_err_long = mean(mean_abs_dev_array(end,:))


%%
m_fig = figure;
hold on
cmap = flipud(brewermap(n_groups,'Spectral'));
colormap(cmap); 
for i = 1:n_groups
    errorbar(t_grid_long,m_mean_array(:,i),m_ste_array(:,i),'Color',cmap(i,:),'Capsize',0,'LineWidth',1.5)
end  

grid on
xlabel('time (burst cycles)');
ylabel('percent deviation (production rate)')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
h = colorbar;
ylabel(h,'average production rate (r)')
ax = gca;
xlim([1e-1 5e3])
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

m_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

m_fig.Renderer='Painters';
set(gca,'xscale','log')
% set(gca,'yscale','log')
ylim([-10 10])
saveas(m_fig,[FigPath,'m_trend_log.png'])
saveas(m_fig,[FigPath,'m_trend_log.pdf'])

% variance
v_fig = figure;
hold on
cmap = flipud(brewermap(n_groups,'Spectral'));
colormap(cmap); 
for i = 1:n_groups
    errorbar(t_grid_long,v_mean_array(:,i),v_ste_array(:,i),'Color',cmap(i,:),'Capsize',0,'LineWidth',1.5)
end  

grid on
xlabel('time (burst cycles)');
ylabel('percent deviation (variance)')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
h = colorbar;
ylabel(h,'average production rate (r)')
ax = gca;
xlim([1e-1 5e3])
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

v_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

v_fig.Renderer='Painters';
set(gca,'xscale','log')
% set(gca,'yscale','log')
ylim([-20 100])
saveas(v_fig,[FigPath,'v_trend_log.png'])
saveas(v_fig,[FigPath,'v_trend_log.pdf'])

