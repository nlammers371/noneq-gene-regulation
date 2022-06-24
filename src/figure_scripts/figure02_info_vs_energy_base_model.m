% Plot results of sharpness parameter sweeps
% clear 
close all
addpath(genpath('../utilities/'))

% %%%%%%%%%%%%%%%%  set relative read and write paths %%%%%%%%%%%%%%%%%%%%
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\';
DataPath = [DropboxFolder  'SweepOutput\sweeps01_info_vs_energy' filesep ];
FigPath = [DropboxFolder '\manuscript\info_vs_energy' filesep];
mkdir(FigPath);

% load dataset with decision time ranges 
load('decision_limit_info.mat','decision_limit_info')

% %%%%%%%%%%%%%%%%  Set plot parameters and constants %%%%%%%%%%%%%%%%%%%%
n_plot = 3e3; % number of points to plot
markerSize = 75; % marker size

% set sweep options
sweep_options = {'n_sim',20,'n_seeds',5,'n_iters_max',50};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%% Plot Phi versus vs IR for simple 4 state gene circuit %%%%%%%%%%%%
% set path 
nStates = 4;
% functionPath = ['../utilities/metricFunctions/n' num2str(nStates) '_OR/'];
functionPath = '../utilities/metricFunctions/symbolic/n004_s01_ns00_g01';
% 'C:\Users\nlamm\projects\noneq-gene-regulation\src\utilities\metricFunctions\symbolic\n004_s01_ns00_g01'
% get metric names 
paramBounds = repmat([-5; 5],1,9);
[~,~,metric_names_sym] = calculateMetricsSym_v2([]);
ir_index = find(strcmp(metric_names_sym,'DecisionRateNorm'));
rate_index = find(strcmp(metric_names_sym,'Production Rate'));
sharpness_index = find(strcmp(metric_names_sym,'Sharpness'));
cycle_time_index = find(strcmp(metric_names_sym,'CycleTime'));
phi_index = find(strcmp(metric_names_sym,'Phi'));

% run symbolic sweep
tic
[sweep_info_sym4, sweep_results_sym4] = ...
                param_sweep_multi_v3([phi_index ir_index],...
                                    functionPath,sweep_options{:},...
                                    'equilibrium_flag',false,...
                                    'TauCycleTime',1,'downsample_output',1,'paramBounds',paramBounds); 
toc

% calculate maximum info rates
metric_array_phi = vertcat(sweep_results_sym4.metric_array);   
rate_array_phi = vertcat(sweep_results_sym4.rate_array);   
phi_vec = metric_array_phi(:,phi_index);
rate_vec = metric_array_phi(:,rate_index);
sharpness_vec = metric_array_phi(:,sharpness_index);
info_vec = metric_array_phi(:,ir_index) * log2(exp(1));% ./ c_factor;
plot_filter = sharpness_vec>=0;% & precision_vec >=0;
plot_indices = randsample(find(plot_filter),n_plot,true,info_vec(plot_filter));

%% make figures
info_vs_flux = figure;
cmap = brewermap(8,'Purples');
hold on

scatter(phi_vec(plot_indices),info_vec(plot_indices),...
      markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',0.75,...
      'MarkerFaceColor',cmap(5,:)); 

% ylim([2 4.25])
xlim([5e-3 3e1])
    
xlabel('energy dissipation rate (k_BT per burst cycle)');
ylabel('information rate (bits per burst cycle)')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

info_vs_flux.InvertHardcopy = 'off';
set(gcf,'color','w');
set(gca,'xscale','log')

saveas(info_vs_flux,[FigPath 'info_vs_flux.png'])
saveas(info_vs_flux,[FigPath 'info_vs_flux.pdf'])                                                                                                                                                                                                                                                                                                                       

nanmax(info_vec)
%% Plot eq and neq decision time bounds for 4 state case

eq_indices = find(sharpness_vec>=0&phi_vec<=1e3&rate_vec>=0.49&rate_vec<=0.51);
[ir_max_eq, ir_eq_i] = sort(info_vec(eq_indices),'descend');
ir_rates_eq = rate_array_phi(eq_indices(ir_eq_i(1:100)),:);
% metric_vec = calculateMetricsSym_v2(ir_rates_eq, sweep_info_sym4);

paramCell = mat2cell(ir_rates_eq,size(ir_rates_eq,1),ones(1,size(ir_rates_eq,2)));
stateProbs = steadyStateVecFunction(paramCell{:});  

% dt_cw_vec = geomean(vertcat(decision_limit_info.cw_ub,decision_limit_info.cw_lb));
dt_ub_vec = decision_limit_info.n_cycles_ub;
% dt_lb_vec = 1e-10*ones(size(dt_cw_vec)); % not worried about lower bound, so exclude from plot
% dt_mean_vec = dt_lb_vec;%mean(vertcat(dt_ub_vec,dt_lb_vec));

% exclude Arabadopsis for now
plot_vec = [1 1 0 0 1]==1;

% close all
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
    phi_filter = 0 <=phi_vec&phi_vec < 1e-8;%phi_filter_vec(i); 
  elseif i < n_lines
    phi_filter = phi_filter_vec(i-1) <= phi_vec&phi_vec < phi_filter_vec(i); 
  else
    phi_filter = true(size(phi_vec));
  end
    
  [neq_max, optimal_neq_index] = max(info_vec.*(phi_filter&plot_filter));

  % extract cycle times
%   cycle_time = metric_array_phi(optimal_neq_index,cycle_time_index);
  optimal_rates = rate_array_flux(optimal_neq_index,:);%*cycle_time/t_cycle; % adjust to proper cycle time

  for e = 1:length(error_vec)
      sweep_info_sym4_temp = sweep_info_sym4;
      sweep_info_sym4_temp.minError = error_vec(e);
      [VNeq, tau_array_neq(e,i),~] = calculateDecisionMetrics(optimal_rates,sweep_info_sym4_temp);    
      tau_array_neq(e,i) = tau_array_neq(e,i);% / cycle_time;
      ir_array_neq(e,i) = VNeq;%*cycle_time;
  end
  tau_array_alt(:,i) = log((1-error_vec)./error_vec) .*(1-2.*error_vec)./VNeq;%*cycle_time;
end
%
close all
decision_time_plot = figure('Position',[100 100 560 420*142.6/152.4]);%;%('Position',[100 100 560*1.15  420]);
% cmap = flipud(brewermap(n_lines,'Spectral'));
cmap2 = brewermap([],'Set2');


topline = repelem(1e5,length(error_vec));
hold on

% plot eq and noneq regimes
fill([error_vec fliplr(error_vec)], [tau_array_neq(:,1)' topline],cmap2(3,:),'FaceAlpha',0.5,'EdgeAlpha',0)
fill([error_vec fliplr(error_vec)], [tau_array_neq(:,end)' fliplr(tau_array_neq(:,1)')],cmap2(2,:),'FaceAlpha',0.5,'EdgeAlpha',0)
cmap = [];

p1 = plot(error_vec,tau_array_neq(:,1),'-','Color',brighten(cmap2(3,:),-.2),'LineWidth',3); 
p2 = plot(error_vec,tau_array_neq(:,end),'Color',brighten(cmap2(2,:),-.2),'LineWidth',3); 
% plot(error_vec,tau_array(:,n_lines)/60,'Color','k','LineWidth',2); 

% legend([p1 p2],'equilibrium limit','nonequilibrium limit','Location','southeast')    
xlabel('error probability (\epsilon)');
ylabel('decision time (burst cycles)');
set(gca,'xdir','reverse')


grid on
set(gca,'FontSize',12)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

ylim([1 1e3])
xlim([0.05 .45])
set(gca,'yscale','log')
% set(gca,'xscale','log')

decision_time_plot.InvertHardcopy = 'off';
set(gcf,'color','w');
saveas(decision_time_plot,[FigPath 'decision_time_plot.png'])
saveas(decision_time_plot,[FigPath 'decision_time_plot.pdf'])

% add decision time limits
for i = find(plot_vec)
    plot(error_vec,repelem(dt_ub_vec(i),length(error_vec)),'-k')
end    

saveas(decision_time_plot,[FigPath 'decision_time_plot_bounds.png'])
saveas(decision_time_plot,[FigPath 'decision_time_plot_bounds.pdf'])

% print results for equilibrium and non-eq decision times at 10%
[~,mi]= min(abs(error_vec-0.32))
dt_eq = tau_array_neq(mi,1)
dt_neq = tau_array_neq(mi,end)

% %% Add decision time bounds for higher-order systems
% DataPath2 = [DropboxFolder  'SweepOutput\sweeps01_info_vs_energy_v3' filesep ];
% 
% % get metric names for numeric sweeps
% [~,~,metric_names_num] = calculateMetricsNumeric_v3([]);
% rate_index_num = find(strcmp(metric_names_num,'ProductionRate'));
% ir_index_num = find(strcmp(metric_names_num,'IR'));
% phi_index_num = find(strcmp(metric_names_num,'Phi'));
% inv_dtime_index_num = find(strcmp(metric_names_num,'InverseDecisionTime'));
% 
% %%%%%%%%%%%%%%%%%%%%5
% % load 5BS results
% % get list of sweep results files with only 1 genera TF reaction
% bs5_sweep_files = dir([DataPath2 'sweep_results_s05*g01*']);
% bs5_info_files = dir([DataPath2 'sweep_info_s05*g01*']);
% 
% % load
% load([DataPath bs5_sweep_files(1).name])
% load([DataPath bs5_info_files(1).name])
% 
% sweep_results5 = sim_results;
% sweep_info5 = sim_info;
% clear sweep_info
% clear sim_info
% 
% %% find best near-equilibrium system
% metric_array5 = sweep_results5.metric_array;
% rate_array5 = sweep_results5.rate_array;
% phi_vec5 = metric_array5(:,phi_index_num);
% dt_vec5 = 1./metric_array5(:,inv_dtime_index_num);
% ir_vec5 = metric_array5(:,ir_index_num);
% rate_vec5 = metric_array5(:,rate_index_num);
% 
% eq_indices = find(phi_vec5>=5e-2&phi_vec5<=1e-1&rate_vec5>=0.02&rate_vec5<=0.98); % NL: this avoids numerical precision issues at very low Phi
% [ir_max, ir_i] = nanmax(ir_vec5(eq_indices));
% 
% error_vec2 = logspace(log10(0.4999),-2,5e2);
% % switch to correct function path
% rmpath(genpath('../utilities/metricFunctions/'));
% fpath = sweep_info5.functionPath;
% src_i = strfind(fpath,'utilities');
% fpath = ['..\' fpath(src_i:end)];
% addpath(genpath(fpath));
% 
% sweep_info5.a1 = 0.98;
% sweep_info5.a0 = 0.02;
% sweep_info5.equilibrium_flag = true;
% 
% ir_rates_eq5 = rate_array5(eq_indices(ir_i),:);
% dt_vec_eq5 = NaN(size(error_vec2));
% 
% for e = 1:length(error_vec2)    
%     sweep_info5.minError = error_vec2(e);
%     metric_vec = calculateMetricsNumeric_v3(ir_rates_eq5, sweep_info5);
%     dt_vec_eq5(e) =  1/metric_vec(inv_dtime_index_num);    
% end
%             
% 
% %%%%%%%%%%%%%%%%%%%%5
% %% load 5LC results
% % get list of sweep results files with only 1 genera TF reaction
% lc5_sweep_files = dir([DataPath2 'sweep_results_s01*g04*']);
% lc5_info_files = dir([DataPath2 'sweep_info_s01*g04*']);
% 
% % load
% load([DataPath lc5_sweep_files(1).name])
% load([DataPath lc5_info_files(1).name])
% 
% sweep_resultsLC5 = sim_results;
% sweep_infoLC5 = sim_info;
% clear sweep_info
% clear sim_info
% 
% %% find best near-equilibrium system
% metric_arrayLC5 = sweep_resultsLC5.metric_array;
% rate_arrayLC5 = sweep_resultsLC5.rate_array;
% phi_vecLC5 = metric_arrayLC5(:,phi_index_num);
% dt_vecLC5 = 1./metric_arrayLC5(:,inv_dtime_index_num);
% ir_vecLC5 = metric_arrayLC5(:,ir_index_num);
% rate_vecLC5 = metric_arrayLC5(:,rate_index_num);
% 
% a_indices = find(rate_vecLC5>=0.1&rate_vecLC5<=0.8); % NL: this avoids numerical precision issues at very low Phi
% [ir_max, ir_i] = sort(ir_vecLC5(a_indices),'descend');
% 
% error_vec2 = logspace(log10(0.4999),-2,5e2);
% % switch to correct function path
% rmpath(genpath('../utilities/metricFunctions/'));
% fpath = sweep_infoLC5.functionPath;
% src_i = strfind(fpath,'utilities');
% fpath = ['..\' fpath(src_i:end)];
% addpath(genpath(fpath));
% 
% sweep_infoLC5.a1 = 0.98;
% sweep_infoLC5.a0 = 0.02;
% sweep_infoLC5.numerical_precision = 5;
% 
% ir_rates_LC5 = rate_arrayLC5(a_indices(ir_i(1)),:);
% dt_vec_LC5 = NaN(size(error_vec2));
% 
% for e = 1:length(error_vec2)    
%     sweep_infoLC5.minError = error_vec2(e);
%     metric_vec = calculateMetricsNumeric_v3(ir_rates_LC5, sweep_infoLC5);
%     dt_vec_LC5(e) =  1/metric_vec(inv_dtime_index_num);    
% end
% %%
% valCellInit = mat2cell( ir_rates_LC5,size( ir_rates_LC5,1),ones(1,size(ir_rates_LC5,2)));    
% 
% % get rate arrays
% Q_num_cs = RSymFun(valCellInit{:});            
% cs_probs_temp = calculate_ss_num(Q_num_cs,5);
%                             
% %% Now overalay results onto simple model fig
% cmap_gra = brewermap(8,'Greys');
% cmap_pu = brewermap(8,'Purples');
% 
% close all
% decision_time_plot = figure('Position',[100 100 560 420*142.6/152.4]);%;%('Position',[100 100 560*1.15  420]);
% 
% hold on
% 
% % plot eq and noneq regimes
% fill([error_vec fliplr(error_vec)], [tau_array_neq(:,1)' topline],cmap2(3,:),'FaceAlpha',0.5,'EdgeAlpha',0)
% fill([error_vec fliplr(error_vec)], [tau_array_neq(:,end)' fliplr(tau_array_neq(:,1)')],cmap2(2,:),'FaceAlpha',0.5,'EdgeAlpha',0)
% cmap = [];
% 
% p1 = plot(error_vec,tau_array_neq(:,1),'-','Color',brighten(cmap2(3,:),-.2),'LineWidth',3); 
% p2 = plot(error_vec,tau_array_neq(:,end),'Color',brighten(cmap2(2,:),-.2),'LineWidth',3); 
% % plot(error_vec,tau_array(:,n_lines)/60,'Color','k','LineWidth',2); 
% 
% % plot higher-order results
% plot(error_vec2, dt_vec_LC5, '-.', 'Color',brighten(cmap_pu(2+4,:),-0.5),'LineWidth',2); 
% plot(error_vec2, dt_vec_eq5, '--', 'Color',brighten(cmap_gra(5,:),-0.5),'LineWidth',2); 
% 
% % legend([p1 p2],'equilibrium limit','nonequilibrium limit','Location','southeast')    
% xlabel('error probability (\epsilon)');
% ylabel('decision time (burst cycles)');
% set(gca,'xdir','reverse')
% 
% 
% grid on
% set(gca,'FontSize',12)
% set(gca,'Color',[228,221,209]/255) 
% 
% ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.XAxis(1).Color = 'k';
% 
% ylim([1 1e3])
% xlim([0.05 .45])
% set(gca,'yscale','log')
% % set(gca,'xscale','log')
% 
% decision_time_plot.InvertHardcopy = 'off';
% set(gcf,'color','w');
% saveas(decision_time_plot,[FigPath 'decision_time_plot_ho.png'])
% saveas(decision_time_plot,[FigPath 'decision_time_plot_ho.pdf'])
% 
% % add decision time limits
% for i = find(plot_vec)
%     plot(error_vec,repelem(dt_ub_vec(i),length(error_vec)),'-k')
% end    
% 
% saveas(decision_time_plot,[FigPath 'decision_time_plot_bounds_ho.png'])
% saveas(decision_time_plot,[FigPath 'decision_time_plot_bounds_ho.pdf'])