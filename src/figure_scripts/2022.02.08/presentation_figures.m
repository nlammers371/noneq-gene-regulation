% Plot results of info vs energy for presentation purposes

clear 
close all
addpath(genpath('../utilities/'))

% %%%%%%%%%%%%%%%%  set relative read and write paths %%%%%%%%%%%%%%%%%%%%
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\';
DataPath = [DropboxFolder  'SweepOutput\sweeps01_info_vs_energy' filesep ];
FigPath = [DropboxFolder '\presentations\UChicago' filesep];
mkdir(FigPath);

% %%%%%%%%%%%%%%%%  Set plot parameters and constants %%%%%%%%%%%%%%%%%%%%
n_plot = 3e3; % number of points to plot
markerAlpha = 0.5; % marker transparency
markerSize = 75; % marker size

% set sweep options
sweep_options = {'n_sim',10,'n_seeds',5,'n_iters_max',50};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%% Plot Phi versus vs IR for simple 4 state gene circuit %%%%%%%%%%%%
% set path 
nStates = 4;
functionPath = ['../utilities/metricFunctions/n' num2str(nStates) '_OR/'];

% get metric names 
[~,~,metric_names_sym] = calculateMetricsSym_v2([]);
ir_index = find(strcmp(metric_names_sym,'DecisionRateNorm'));
decision_time_index = find(strcmp(metric_names_sym,'DecisionTimeNorm'));
sharpness_index = find(strcmp(metric_names_sym,'Sharpness'));
cycle_time_index = find(strcmp(metric_names_sym,'CycleTime'));
phi_index = find(strcmp(metric_names_sym,'Phi'));

% run symbolic sweep
tic
[sweep_info_sym4, sweep_results_sym4] = ...
                param_sweep_multi_v3([phi_index ir_index],...
                                    functionPath,sweep_options{:},...
                                    'half_max_flag',false,'equilibrium_flag',false,...
                                    'TauCycleTime',1,'downsample_output',1); 
toc
%% calculate maximum info rates
metric_array_phi = vertcat(sweep_results_sym4.metric_array);   
phi_vec = metric_array_phi(:,phi_index);
sharpness_vec = metric_array_phi(:,sharpness_index);
info_vec = metric_array_phi(:,ir_index) * log2(exp(1));% ./ c_factor;
plot_filter = sharpness_vec>=0;% & precision_vec >=0;
plot_indices = randsample(find(plot_filter),n_plot,true,info_vec(plot_filter));

t_cycle = 5;

% make figures
info_vs_flux = figure;
cmap = brewermap(8,'Purples');
hold on

scatter(phi_vec(plot_indices),info_vec(plot_indices)/t_cycle,...
      markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',0.75,...
      'MarkerFaceColor',cmap(5,:)); 

% ylim([2 4.25])
xlim([5e-3 3e1])
    
xlabel('energy dissipation rate (k_BT per minute)');
ylabel('information rate (bits per minute)')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
set(gca,'xtick',[1e-2 1e-1 1e0 1e1])

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

info_vs_flux.InvertHardcopy = 'off';
set(gcf,'color','w');
set(gca,'xscale','log')

saveas(info_vs_flux,[FigPath 'info_vs_flux.png'])
saveas(info_vs_flux,[FigPath 'info_vs_flux.pdf'])                                                                                                                                                                                                                                                                                                                       

%% Plot eq and neq decision time bounds for 4 state case

% close all
% t_cycle = 10*60; % set the cycle time
n_lines = 10;
error_vec = logspace(log10(0.4999),-2,1e3);

% generate vector of energy levels
phi_filter_vec = logspace(-3,2,n_lines);
rate_array_flux = vertcat(sweep_results_sym4.rate_array); 

% initialize arrays to store rates
rate_array = NaN(n_lines,size(rate_array_flux,2));
tau_array_neq = NaN(length(error_vec),n_lines);
tau_array_alt = NaN(length(error_vec),n_lines);
ir_array_neq = NaN(length(error_vec),n_lines);

% iterate
for i = 1:n_lines
  if i == 1 % for equilibrium
    phi_filter = 0 <=phi_vec&phi_vec < phi_filter_vec(i); 
  else
    phi_filter = phi_filter_vec(i-1) <= phi_vec&phi_vec < phi_filter_vec(i); 
  end
    
  [neq_max, optimal_neq_index] = max(info_vec.*(phi_filter&plot_filter));

  % extract cycle times
  cycle_time = metric_array_phi(optimal_neq_index,cycle_time_index);
  optimal_rates = rate_array_flux(optimal_neq_index,:);%*cycle_time/t_cycle; % adjust to proper cycle time

  for e = 1:length(error_vec)
      sweep_info_sym4_temp = sweep_info_sym4;
      sweep_info_sym4_temp.minError = error_vec(e);
      [VNeq,tau_array_neq(e,i),~] = calculateDecisionMetrics(optimal_rates,sweep_info_sym4_temp);    
      tau_array_neq(e,i) = tau_array_neq(e,i) / cycle_time;
      ir_array_neq(e,i) = VNeq*cycle_time;
  end
  tau_array_alt(:,i) = log((1-error_vec)./error_vec) .*(1-2.*error_vec)./VNeq*cycle_time;
end
  
close all
decision_time_plot = figure;%('Position',[100 100 560*1.15  420]);
% cmap = flipud(brewermap(n_lines,'Spectral'));
cmap2 = brewermap([],'Set2');


topline = repelem(1e5,length(error_vec));
hold on

% plot eq and noneq regimes
fill([error_vec fliplr(error_vec)], [tau_array_neq(:,1)' topline]*t_cycle,cmap2(3,:),'FaceAlpha',0.5,'EdgeAlpha',0)
fill([error_vec fliplr(error_vec)], [tau_array_neq(:,end)' fliplr(tau_array_neq(:,1)')]*t_cycle,cmap2(2,:),'FaceAlpha',0.5,'EdgeAlpha',0)
cmap = [];
% plot lines
% for i = 2:n_lines
% %     plot(error_vec,tau_array(:,i)/60,'Color',[cmap(i,:) 1],'LineWidth',1); 
%     cmap(i-1,:) = brighten(cmap2(2,:),.9-.16*i);  
%     fill([error_vec fliplr(error_vec)], [tau_array(:,i-1)'/60 fliplr(tau_array(:,i)')/60],cmap(i-1,:),'FaceAlpha',0.6,'EdgeAlpha',.2)
% end
% colormap(vertcat(cmap2(3,:),cmap))
p1 = plot(error_vec,tau_array_neq(:,1)*t_cycle,'-','Color',brighten(cmap2(3,:),-.2),'LineWidth',3); 
p2 = plot(error_vec,tau_array_neq(:,end)*t_cycle,'Color',brighten(cmap2(2,:),-.2),'LineWidth',3); 
% plot(error_vec,tau_array(:,n_lines)/60,'Color','k','LineWidth',2); 

legend([p1 p2],'equilibrium limit','nonequilibrium limit','Location','southeast')    
xlabel('error probability (\epsilon)');
ylabel('decision time (minutes)');
set(gca,'xdir','reverse')

% h = colorbar;
% h.Ticks = [.05 round(linspace(1.5,9.5,4)/10,1)];
% h.TickLabels = {'0','10^{-2}','10^{-1}','10^{0}','10^{1}'};
% h.Color = 'k';
% ylabel(h,'\Phi (k_BT per cycle)')
% set(gca,'xtick',[0.05 0.15 
grid on
set(gca,'FontSize',12)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

ylim([5 3e3]*t_cycle)
xlim([0.05 .45])
set(gca,'yscale','log')
% set(gca,'xscale','log')

decision_time_plot.InvertHardcopy = 'off';
set(gcf,'color','w');
saveas(decision_time_plot,[FigPath 'decision_time_plot.png'])
saveas(decision_time_plot,[FigPath 'decision_time_plot.pdf'])

% print results for equilibrium and non-eq decision times at 10%
[~,mi]= min(abs(error_vec-0.32))
dt_eq = tau_array_neq(mi,1)
dt_neq = tau_array_neq(mi,end)

%% Plot DT vs Energy

% generate vector of energy levels
phi_bins = logspace(-3,log10(35),51);
rate_array_flux = vertcat(sweep_results_sym4.rate_array); 
minError = 0.32;

% initialize arrays to store rates
dt_bound_vec = NaN(1,length(phi_bins)-1);
ir_bound_vec = NaN(1,length(phi_bins)-1);

% iterate
for i = 1:length(phi_bins)-1
 
  % identify valid systems
  phi_filter = phi_bins(i)<=phi_vec & phi_vec<phi_bins(i+1);  
    
  % find best gene circuit
  [neq_max, optimal_neq_index] = nanmax(info_vec.*(phi_filter));

  % extract cycle times  
  optimal_rates = rate_array_flux(optimal_neq_index,:);
  sweep_info_sym4_temp = sweep_info_sym4;
  sweep_info_sym4_temp.minError = minError;
  [ir_bound_vec(i),dt_bound_vec(i),~] = calculateDecisionMetrics(optimal_rates,sweep_info_sym4_temp);    
    
end


decision_time_plot = figure;
cmap = brewermap(8,'Purples');

phi_axis = phi_bins(1:end-1) + diff(phi_bins);
topline = repelem(1e5,length(phi_axis));

hold on

% plot eq and noneq regimes
fill([phi_axis fliplr(phi_axis)], [imgaussfilt(dt_bound_vec,1) topline]*t_cycle, cmap(5,:),'FaceAlpha',0.5,'EdgeAlpha',0)


p1 = plot(phi_axis,imgaussfilt(dt_bound_vec,1)*t_cycle,'-','Color','k','LineWidth',3); 

xlabel('energy dissipation rate (k_BT per minute)');
ylabel('decision time (minutes)');

grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
set(gca,'xtick',[1e-2 1e-1 1e0 1e1])
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

ylim([0 600])
xlim([5e-3 3e1])

% set(gca,'yscale','log')
set(gca,'xscale','log')

decision_time_plot.InvertHardcopy = 'off';
set(gcf,'color','w');
saveas(decision_time_plot,[FigPath 'decision_time_vs_phi.png'])
saveas(decision_time_plot,[FigPath 'decision_time_vs_phi.pdf'])


%% make figures
info_vs_flux_bound = figure;
cmap = brewermap(8,'Purples');
hold on

scatter(phi_vec(plot_indices),info_vec(plot_indices)/t_cycle,...
      markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',0.75,...
      'MarkerFaceColor',cmap(5,:)); 

plot(phi_axis,ir_bound_vec/t_cycle*log2(exp(1)),'-k','LineWidth',3)    
% ylim([2 4.25])
xlim([5e-3 3e1])
    
xlabel('energy dissipation rate (k_BT per minute)');
ylabel('information rate (bits per minute)')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
set(gca,'xtick',[1e-2 1e-1 1e0 1e1])

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

info_vs_flux_bound.InvertHardcopy = 'off';
set(gcf,'color','w');
set(gca,'xscale','log')

saveas(info_vs_flux_bound,[FigPath 'info_vs_flux_bound.png'])
saveas(info_vs_flux_bound,[FigPath 'info_vs_flux_bound.pdf'])   