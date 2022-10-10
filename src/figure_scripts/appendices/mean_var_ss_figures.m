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
load([ReadPath 'Gaussian_noise_struct.mat'])

%% Calculate percent absolute deviation from ground truth
var_dev_vec = NaN(length(Gaussian_noise_struct(1).r_mean_vec),length(Gaussian_noise_struct));
mean_dev_vec = NaN(length(Gaussian_noise_struct(1).r_mean_vec),length(Gaussian_noise_struct));

var_abs_dev_vec = NaN(length(Gaussian_noise_struct(1).r_mean_vec),length(Gaussian_noise_struct));
mean_abs_dev_vec = NaN(length(Gaussian_noise_struct(1).r_mean_vec),length(Gaussian_noise_struct));

for i = 1:length(Gaussian_noise_struct)
    rt = Gaussian_noise_struct(i).r_predicted;
    vt = Gaussian_noise_struct(i).var_predicted;

    mean_dev_vec(:,i) = (rt-Gaussian_noise_struct(i).r_mean_vec)/rt;
    var_dev_vec(:,i) = (vt-Gaussian_noise_struct(i).r_var_vec)/vt;

    mean_abs_dev_vec(:,i) = abs(rt-Gaussian_noise_struct(i).r_mean_vec)/rt;
    var_abs_dev_vec(:,i) = abs(vt-Gaussian_noise_struct(i).r_var_vec)/vt;

end

%% now look at p values for gene circuits, grouped by mean rate
n_groups = 10;
r_bins = linspace(0,1,n_groups+1);
r_mean_vec = [Gaussian_noise_struct.r_mean];
r_group_vec = discretize(r_mean_vec,r_bins);
% p_array_sim = p_array';%vertcat(Gaussian_noise_struct.p_vec);

% use bootstrapping to estimate standard error in convergence across
% different groups
t_grid_long = Gaussian_noise_struct(1).time_grid;
nBoots = 100;
ma_boot_array = NaN(length(t_grid_long),n_groups,nBoots);
m_boot_array = NaN(length(t_grid_long),n_groups,nBoots);
va_boot_array = NaN(length(t_grid_long),n_groups,nBoots);
v_boot_array = NaN(length(t_grid_long),n_groups,nBoots);
for n = 1:n_groups
    option_vec = find(r_group_vec==n);
    for b = 1:nBoots
        boot_indices = randsample(option_vec,length(option_vec),true);
       
        ma_boot_array(:,n,b) = mean(mean_abs_dev_vec(:,boot_indices),2);
        m_boot_array(:,n,b) = mean(mean_dev_vec(:,boot_indices),2);
        va_boot_array(:,n,b) = mean(var_abs_dev_vec(:,boot_indices),2);
        v_boot_array(:,n,b) = mean(var_dev_vec(:,boot_indices),2);
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
ylim([0 30])
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
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

va_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

va_fig.Renderer='Painters';
set(gca,'xscale','log')
% set(gca,'yscale','log')
ylim([0 100])
saveas(va_fig,[FigPath,'va_trend_log.png'])
saveas(va_fig,[FigPath,'va_trend_log.pdf'])


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
va_fig = figure;
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

ma_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

ma_fig.Renderer='Painters';
set(gca,'xscale','log')
% set(gca,'yscale','log')
ylim([-20 100])
saveas(ma_fig,[FigPath,'v_trend_log.png'])
saveas(ma_fig,[FigPath,'v_trend_log.pdf'])
