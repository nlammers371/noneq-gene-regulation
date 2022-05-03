% script to plot simulated bursting trends
clear 
close all
addpath(genpath('../utilities/'))
                        
% define read path
readPath = ['../../out/illustrative_bursting_simulations' filesep];

% define save path
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\';
DataPath = [DropboxFolder  'SweepOutput\sweeps01_info_vs_energy' filesep ];
figPath = [DropboxFolder '\manuscript\illustrative_bursting_simulations' filesep];
mkdir(figPath);

% load data 
load([readPath 'sim_struct.mat']);
time_grid = sim_struct(1).time_grid;

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Perform some basic calculations

%%%%%%%%%%
%% NEQ
x_max = 50; % burst cycles
last_ind = find(time_grid==x_max);

% calculations for c0
mRNA_c0_mean_neq = mean(sim_struct(1).total_mRNA_array,2)';
% find trace that is most representative of mean
[~,mi_c0_neq] = min(sum(abs(sim_struct(1).total_mRNA_array(1:last_ind,:)-mRNA_c0_mean_neq(1:last_ind)'),1),[],2);
mRNA_c0_rep_neq = sim_struct(1).total_mRNA_array(:,mi_c0_neq)';
mRNA_c0_ste_neq = std(sim_struct(1).total_mRNA_array,[],2)';
mRNA_c0_ub_neq = mRNA_c0_rep_neq + mRNA_c0_ste_neq;
mRNA_c0_lb_neq = mRNA_c0_rep_neq - mRNA_c0_ste_neq;

% calculations for c1
mRNA_c1_mean_neq = mean(sim_struct(2).total_mRNA_array,2)';
[~,mi_c1_neq] = min(sum(abs(sim_struct(2).total_mRNA_array(1:last_ind,:)-mRNA_c1_mean_neq(1:last_ind)'),1),[],2);
mRNA_c1_rep_neq = sim_struct(2).total_mRNA_array(:,mi_c1_neq)';
mRNA_c1_ste_neq = std(sim_struct(2).total_mRNA_array,[],2)';
mRNA_c1_ub_neq = mRNA_c1_rep_neq + mRNA_c1_ste_neq;
mRNA_c1_lb_neq = mRNA_c1_rep_neq - mRNA_c1_ste_neq;

%%%%%%%%%%
%%% EQ
% calculations for c0
mRNA_c0_mean_eq = mean(sim_struct(3).total_mRNA_array,2)';
% find trace that is most representative of mean
[~,mi_c0_eq] = min(sum(abs(sim_struct(3).total_mRNA_array(1:last_ind,:)-mRNA_c0_mean_eq(1:last_ind)'),1),[],2);
mRNA_c0_rep_eq = sim_struct(3).total_mRNA_array(:,mi_c0_eq)';
mRNA_c0_ste_eq = std(sim_struct(3).total_mRNA_array,[],2)';
mRNA_c0_ub_eq = mRNA_c0_rep_eq + mRNA_c0_ste_eq;
mRNA_c0_lb_eq = mRNA_c0_rep_eq - mRNA_c0_ste_eq;

% calculations for c1
mRNA_c1_mean_eq = mean(sim_struct(4).total_mRNA_array,2)';
[~,mi_c1_eq] = min(sum(abs(sim_struct(4).total_mRNA_array(1:last_ind,:)-mRNA_c1_mean_eq(1:last_ind)'),1),[],2);
mRNA_c1_rep_eq = sim_struct(4).total_mRNA_array(:,mi_c1_eq)';
mRNA_c1_ste_eq = std(sim_struct(4).total_mRNA_array,[],2)';
mRNA_c1_ub_eq = mRNA_c1_rep_eq + mRNA_c1_ste_eq;
mRNA_c1_lb_eq = mRNA_c1_rep_eq - mRNA_c1_ste_eq;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot simulated burst trajectories
close all
burst_x = 25;
b_factor = 0.25;
% downsample for 
% add noise for illustrative purposes
p_path_c0_neq = sim_struct(1).promoter_state_array(:,mi_c0_neq)';
n0_vec = imgaussfilt(normrnd(0,0.05,1,length(p_path_c0_neq)),1);

p_path_c1_neq = sim_struct(2).promoter_state_array(:,mi_c1_neq)';
n1_vec = imgaussfilt(normrnd(0,0.05,1,length(p_path_c1_neq)),1);

%%%%%%% C0 NEQ
neq_burst_fig_c0 = figure('Position',[100 100 512 128]);
cmap = brewermap([],'Set2');
hold on
% plot sample trajectory 
plot(time_grid,p_path_c0_neq+n0_vec, 'Color', brighten(cmap(2,:),b_factor), 'LineWidth',1);

% formatting
xlabel('time');
ylabel('transcription rate')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
xlim([0 burst_x])
ylim([-0.1 1.1])
grid on
neq_burst_fig_c0.InvertHardcopy = 'off';
set(gcf,'color','w');
set(gca,'YTick',[0 1])
box on
saveas(neq_burst_fig_c0,[figPath 'burst_trend_neq_c0.png'])
saveas(neq_burst_fig_c0,[figPath 'burst_trend_neq_c0.pdf'])

%%%%%%% C1 NEQ
neq_burst_fig_c1 = figure('Position',[100 100 512 128]);

hold on
% plot sample trajectory 
plot(time_grid,p_path_c1_neq+n1_vec, 'Color', brighten(cmap(2,:),-b_factor), 'LineWidth',1);

% formatting
xlabel('time');
ylabel('transcription rate')

set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
xlim([0 burst_x])
ylim([-0.1 1.1])

neq_burst_fig_c1.InvertHardcopy = 'off';
set(gcf,'color','w');
set(gca,'YTick',[0 1])
box on
saveas(neq_burst_fig_c1,[figPath 'burst_trend_neq_c1.png'])
saveas(neq_burst_fig_c1,[figPath 'burst_trend_neq_c1.pdf'])

%%%%%%% C0 EQ

% add noise for illustrative purposes
p_path_c0_eq = sim_struct(3).promoter_state_array(:,mi_c0_eq)';
n0_vec = imgaussfilt(normrnd(0,0.05,1,length(p_path_c0_eq)),1);

p_path_c1_eq = sim_struct(4).promoter_state_array(:,mi_c1_eq)';
n1_vec = imgaussfilt(normrnd(0,0.05,1,length(p_path_c1_eq)),1);

eq_burst_fig_c0 = figure('Position',[100 100 512 128]);
hold on
% plot sample trajectory 
plot(time_grid,p_path_c0_eq+n0_vec, 'Color', brighten(cmap(3,:),b_factor), 'LineWidth',1);

% formatting
xlabel('time');
ylabel('transcription rate')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
xlim([0 burst_x])
ylim([-0.1 1.1])

eq_burst_fig_c0.InvertHardcopy = 'off';
set(gcf,'color','w');
set(gca,'YTick',[0 1])
box on
saveas(eq_burst_fig_c0,[figPath 'burst_trend_eq_c0.png'])
saveas(eq_burst_fig_c0,[figPath 'burst_trend_eq_c0.pdf'])

%%%%%%% C1 NEQ
eq_burst_fig_c1 = figure('Position',[100 100 512 128]);

hold on
% plot sample trajectory 
plot(time_grid,p_path_c1_eq+n1_vec, 'Color', brighten(cmap(3,:),-b_factor), 'LineWidth',1);

% formatting
xlabel('time');
ylabel('transcription rate')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
xlim([0 burst_x])
ylim([-0.1 1.1])

eq_burst_fig_c1.InvertHardcopy = 'off';
set(gcf,'color','w');
set(gca,'YTick',[0 1])
box on
saveas(eq_burst_fig_c1,[figPath 'burst_trend_eq_c1.png'])
saveas(eq_burst_fig_c1,[figPath 'burst_trend_eq_c1.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make accumulated mRNA figures
neq_mRNA_fig = figure('Position',[100 100 512 200]);
hold on
x_max = 75;
% plot uncertainty cones
fill([time_grid fliplr(time_grid)],[mRNA_c0_lb_neq fliplr(mRNA_c0_ub_neq)],brighten(cmap(2,:),b_factor),'EdgeAlpha',0,'FaceAlpha',0.5)
fill([time_grid fliplr(time_grid)],[mRNA_c1_lb_neq fliplr(mRNA_c1_ub_neq)],brighten(cmap(2,:),-b_factor),'EdgeAlpha',0,'FaceAlpha',0.5)

% plot sample trajectory 
p0 = plot(time_grid, mRNA_c0_rep_neq, 'Color', brighten(cmap(2,:),b_factor), 'LineWidth',2);
p1 = plot(time_grid, mRNA_c1_rep_neq, 'Color', brighten(cmap(2,:),-b_factor), 'LineWidth',2);

% formatting
xlabel('time');
ylabel('accumulated mRNA')

set(gca,'FontSize',14,'xtick',0:25:75)
% legend([p0 p1],'c_0','c_1','Location','northwest')
set(gca,'Color',[228,221,209]/255) 
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
xlim([0 x_max])
ylim([0 50])

% grid on
neq_mRNA_fig.InvertHardcopy = 'off';
neq_mRNA_fig.Renderer = 'Painters';
set(gcf,'color','w');
box on 
saveas(neq_mRNA_fig,[figPath 'mRNA_trend_neq.png'])
saveas(neq_mRNA_fig,[figPath 'mRNA_trend_neq.pdf'])



eq_mRNA_fig = figure('Position',[100 100 512 200]);
hold on

% plot uncertainty cones
fill([time_grid fliplr(time_grid)],[mRNA_c0_lb_eq fliplr(mRNA_c0_ub_eq)],brighten(cmap(3,:),b_factor),'EdgeAlpha',0,'FaceAlpha',0.5)
fill([time_grid fliplr(time_grid)],[mRNA_c1_lb_eq fliplr(mRNA_c1_ub_eq)],brighten(cmap(3,:),-b_factor),'EdgeAlpha',0,'FaceAlpha',0.5)

% plot sample trajectory 
p0 = plot(time_grid, mRNA_c0_rep_eq, 'Color', brighten(cmap(3,:),b_factor), 'LineWidth',2);
p1 = plot(time_grid, mRNA_c1_rep_eq, 'Color', brighten(cmap(3,:),-b_factor), 'LineWidth',2);

% formatting
xlabel('time');
ylabel('accumulated mRNA')

set(gca,'FontSize',14,'xtick',0:25:75)
% legend([p0 p1],'c_0','c_1','Location','northwest')
set(gca,'Color',[228,221,209]/255) 
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
xlim([0 x_max])
ylim([0 50])
% grid on
eq_mRNA_fig.Renderer = 'Painters';
eq_mRNA_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
box on
saveas(eq_mRNA_fig,[figPath 'mRNA_trend_eq.png'])
saveas(eq_mRNA_fig,[figPath 'mRNA_trend_eq.pdf'])