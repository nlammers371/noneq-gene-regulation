% This script to make figures examining validity of Gaussian noise
% approximation

clear
close all
addpath(genpath('../utilities'))
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\';
ReadPath = [DropboxFolder  'SweepOutput\appendices' filesep ];
FigPath = [DropboxFolder 'manuscript\appendices\gaussian_approximation\'];
mkdir(FigPath)
% load simulation data
load([ReadPath 'Gaussian_noise_struct.mat'])

%% make figures
close all
% plot accumulated mRNA over time
sim_ind = 11;
mRNA_array = Gaussian_noise_struct(sim_ind).mRNA_array;
time_grid = Gaussian_noise_struct(sim_ind).time_grid;

mRNA_fig = figure;
hold on
cmap = brewermap(size(mRNA_array,2),'Greens');
rng(476);
for i = randsample(1:size(mRNA_array,2),size(mRNA_array,2),false)
    plot(time_grid,mRNA_array(:,i)./(1+time_grid'),'Color',cmap(i,:))
end  

grid on
xlabel('time (burst cycles)');
ylabel('average transcription rate (r)')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

mRNA_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
mRNA_fig.Renderer='Painters';

saveas(mRNA_fig,[FigPath,'mRNA_plot.png'])
saveas(mRNA_fig,[FigPath,'mRNA_plot.pdf'])

%% Plot illustrative histograms from 3 different time points
plot_times = [10 50 5000];
plot_indices = find(ismember(time_grid,plot_times));
cmap = brewermap(length(plot_indices)+2,'Greens');
% hist_bins = linspace(
for p = 1:length(plot_indices)
  
    hist_fig = figure;
    histogram(mRNA_array(plot_indices(p),:)/(1+time_grid(plot_indices(p))),'Normalization','probability','FaceColor',cmap(p+1,:))
    xlabel('transcription rate');
    ylabel('probability')
    grid on
    set(gca,'FontSize',14)
    set(gca,'Color',[228,221,209]/255) 

    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';
    
    hist_fig.InvertHardcopy = 'off';
    set(gcf,'color','w');
    xlim([0 1])
    ylim([0 0.4])
    saveas(hist_fig,[FigPath,'mRNA_hist_t' sprintf('%04d',plot_times(p)) '.png'])
    saveas(hist_fig,[FigPath,'mRNA_hist_t' sprintf('%04d',plot_times(p)) '.pdf'])
end  

%% now look at p values for gene circuits, grouped by mean rate
n_groups = 10;
r_bins = linspace(0,1,n_groups+1);
r_mean_vec = [Gaussian_noise_struct.r_mean];
r_group_vec = discretize(r_mean_vec,r_bins);
p_array_sim = vertcat(Gaussian_noise_struct.p_vec);

% use bootstrapping to estimate standard error in convergence across
% different groups
t_grid_long = 0:max(time_grid);
nBoots = 100;
p_boot_array = NaN(length(t_grid_long),n_groups,nBoots);
for n = 1:n_groups
    option_vec = find(r_group_vec==n);
    for b = 1:nBoots
        boot_indices = randsample(option_vec,length(option_vec),true);
        p_array_boot = p_array_sim(boot_indices,:);
        p_boot_array(:,n,b) = nanmean(p_array_boot,1)';
    end
end

p_mean_array = nanmean(p_boot_array,3);
p_ste_array = nanstd(p_boot_array,[],3);

%%
close all
p_fig = figure;
hold on
cmap = flipud(brewermap(n_groups,'Spectral'));
colormap(cmap); 
for i = 1:n_groups
    errorbar(t_grid_long+1,p_mean_array(:,i),p_ste_array(:,i),'Color',cmap(i,:),'Capsize',0,'LineWidth',1.5)
end  

grid on
xlabel('time (burst cycles)');
ylabel('p-value (kstest)')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
h = colorbar;
ylabel(h,'average production rate (r)')
ax = gca;

ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

p_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

% saveas(p_fig,[FigPath,'p_trend.png'])
% saveas(p_fig,[FigPath,'p_trend.pdf'])
p_fig.Renderer='Painters';
set(gca,'xscale','log')

saveas(p_fig,[FigPath,'p_trend_log.png'])
saveas(p_fig,[FigPath,'p_trend_log.pdf'])

