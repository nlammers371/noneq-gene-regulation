% script to track information rate as a function of cw

clear 
close all
addpath(genpath('../utilities/'))

% set basic parameters
nStates = 6;
paramBounds = repmat([-5 ; 5],1,11); % constrain transition rate magnitude
[~,~,metric_names] = calculateMetricsSym_v2([]);

% specify function path
% functionPath = ['../utilities/metricFunctions/n' num2str(nStates) '_OR/'];
functionPath = '../utilities/metricFunctions/symbolic/n006_s01_ns01_g01';

% make sure we're linked to the appropriate function subfolder% make sure we're linked to the appropriate function subfolder
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));

% load dataset with decision time ranges 
load('decision_limit_info.mat','decision_limit_info')

% define save path
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\manuscript\';
FigPath = [DropboxFolder 'info_vs_cw' filesep];
mkdir(FigPath);         

% get index of useful metrics
ir_index = find(strcmp(metric_names,'DecisionRateNorm'));
rate_index = find(strcmp(metric_names,'Production Rate'));
phi_index = find(strcmp(metric_names,'Phi'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
cw_index = find(strcmp(metric_names,'CW'));
cycle_time_index = find(strcmp(metric_names,'CycleTime'));
tau_index = find(strcmp(metric_names,'DecisionTimeNorm'));

% set sim options
sweep_options = {'n_seeds',10,'n_iters_max',50,'n_sim',20,'nStates',nStates};

% calculate sensitivity bound
alpha_factor = 100; % note that I employ a lower alpha factor here so permit exploration of higher cw concentrations. Results only depend on relative alpha scale
% seq = 1/4;

% specify plot options 
markerAlpha = 0.5; % marker transparency
markerSize = 75; % marker size
TauCycleTime = 1; % do this for circuits with a 5 minute cycle time

bit_factor = log2(exp(1));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% info vs cw %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
[sim_info_neq, sim_struct_neq] = param_sweep_multi_v3([ir_index cw_index],functionPath,sweep_options{:},...
                                          'TauCycleTime',TauCycleTime,...
                                          'equilibrium_flag',false,'specFactor',alpha_factor,'paramBounds',paramBounds);

[sim_info_eq, sim_struct_eq] = param_sweep_multi_v3([ir_index cw_index],functionPath,sweep_options{:},...
                                          'TauCycleTime',TauCycleTime,...
                                          'equilibrium_flag',true,'specFactor',alpha_factor,'paramBounds',paramBounds);                                        
toc     

%% pull result arrays
n_plot = 3e3;
metric_array_neq = vertcat(sim_struct_neq.metric_array);
metric_array_eq = vertcat(sim_struct_eq.metric_array);

% extract key vectors
sharpness_vec_neq = metric_array_neq(:,sharpness_index);
sharpness_vec_eq = metric_array_eq(:,sharpness_index);

ir_vec_neq = metric_array_neq(:,ir_index)*bit_factor;
ir_vec_eq = metric_array_eq(:,ir_index)*bit_factor;

cw_vec_neq = metric_array_neq(:,cw_index);
cw_vec_eq = metric_array_eq(:,cw_index);

neq_indices = find(sharpness_vec_neq > 0);
eq_indices = find(sharpness_vec_eq > 0);


% Make plots
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First plot just eq vs. neq 

% generate vectors to plot (neq)
plot_indices_neq = randsample(neq_indices,min([n_plot,length(neq_indices)]),false);

% generate vectors to plot (eq)
plot_indices_eq = randsample(eq_indices,min([n_plot,length(eq_indices)]),false);

% calculate upper IR bounds for EQ and NEQ
cw_axis = logspace(0,6);
neq_bound_vec = NaN(1,length(cw_axis)-1);
eq_bound_vec = NaN(1,length(cw_axis)-1);

for c = 1:length(cw_axis)-1
    cw_neq_filter = 10.^cw_vec_neq(neq_indices)>=cw_axis(c) & 10.^cw_vec_neq(neq_indices)<cw_axis(c+1);
    if any(cw_neq_filter)
      neq_bound_vec(c) = nanmax(ir_vec_neq(neq_indices(cw_neq_filter)));
    end
    cw_eq_filter = 10.^cw_vec_eq(eq_indices)>=cw_axis(c) & 10.^cw_vec_eq(eq_indices)<cw_axis(c+1);
    if any(cw_eq_filter)
      eq_bound_vec(c) = nanmax(ir_vec_eq(eq_indices(cw_eq_filter)));
    end
end
% cw_axis = cw_axis(1:end-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make figure
ir_vs_cw_eq_fig = figure('Position',[100 100 560 420]);
hold on
cmap = brewermap([],'Set2');

sneq = scatter(10.^cw_vec_neq(plot_indices_neq), ir_vec_neq(plot_indices_neq),...
      markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',0.25, 'MarkerFaceColor',cmap(2,:));
    
seq = scatter(10.^cw_vec_eq(plot_indices_eq), ir_vec_eq(plot_indices_eq),...
      markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',0.25, 'MarkerFaceColor',cmap(3,:));    
    
% plot(10.^cw_axis,imgaussfilt(neq_bound_vec),'--','Color',brighten(cmap(2,:),-.5),'LineWidth',3)
% plot(10.^cw_axis,imgaussfilt(eq_bound_vec),'--','Color',brighten(cmap(3,:),-.5),'LineWidth',3)

set(gca,'xscale','log') 
set(gca,'yscale','log') 
xlabel('non-cognate factor concentration (w/c)');
ylabel('information rate (bits per cycle)')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
set(gca,'xtick',[1 10 100 10^3 10^4 10^5]);%,'xticklabels',{'\alpha^0','\alpha^1','\alpha^{2}','\alpha^3','\alpha^4'})
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
ylim([1e-8 .1])
xlim([1 1e5])
grid on
ir_vs_cw_eq_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(ir_vs_cw_eq_fig,[FigPath 'IR_vs_cw.png'])
saveas(ir_vs_cw_eq_fig,[FigPath 'IR_vs_cw.pdf'])

% plot reference curves
ir_eq = 0.0035;
occ_vec_eq = 100./cw_axis;
% occ_vec_neq = 100^2./cw_axis;

% plot(cw_axis,(occ_vec_eq./(occ_vec_eq+1)).*ir_eq,'-','LineWidth',2,'Color',brighten(cmap(2,:),-0.65))
plot(cw_axis,(occ_vec_eq./(occ_vec_eq+1)).^2.*ir_eq*4,'-.','LineWidth',2,'Color',brighten(cmap(2,:),-0.65))
plot(cw_axis,(occ_vec_eq./(occ_vec_eq+1)).^2.*ir_eq,'-.','LineWidth',2,'Color',brighten(cmap(3,:),-0.65))

% saveas(ir_vs_cw_eq_fig,[FigPath 'IR_vs_cw_with_bounds.png'])
% saveas(ir_vs_cw_eq_fig,[FigPath 'IR_vs_cw_with_bounds.pdf'])



%% Plot illustrative occupancy curve for single binding site at equilibrium
occ_pd_eq = alpha_factor ./ (alpha_factor + cw_axis);

occ_fig = figure('Position',[100 100 560 420]);
hold on

plot(cw_axis, occ_pd_eq, '-k','LineWidth',3)  
    
set(gca,'xscale','log') 

xlabel('non-cognate factor concentration (w/c)');
ylabel('cognate factor occupancy (p_c)')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
set(gca,'xtick',[1 10 10^2 10^3 10^4 10^5]);%,'xticklabels',{'\alpha^0','\alpha^1','\alpha^{2}','\alpha^3','\alpha^4'})
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
ylim([0 1.05])
xlim([1 1e5])
grid on
occ_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(occ_fig,[FigPath 'bs_occupancy.png'])
saveas(occ_fig,[FigPath 'bs_occupancy.pdf'])

%% calculate relative non-eq gain as a function of CW
close all

gain_fig = figure;
peq = plot(cw_axis(1:end-1),imgaussfilt(neq_bound_vec./eq_bound_vec,2),'Color','k','LineWidth',3);
set(gca,'xscale','log') 
set(gca,'yscale','log') 
xlabel('non-cognate factor concentration (w/c)');
ylabel('information gain')
box on
set(gca,'LineWidth',2)
set(gca,'FontSize',20)
% set(gca,'Color',[228,221,209]/255) 
set(gca,'xtick',[1  10^2  10^4 ]);%,'xticklabels',{'\alpha^0','\alpha^1','\alpha^{2}','\alpha^3','\alpha^4'})
% set(gca,'ytick',[1  alpha_factor^1  alpha_factor^2],'yticklabels',{'\alpha^0','\alpha^{1}','\alpha^2'})
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
% ylim([10^-6.7 1])
xlim([1 1e5])
% grid on
gain_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(gain_fig,[FigPath 'neq_ir_gain.png'])
saveas(gain_fig,[FigPath 'neq_ir_gain.pdf'])

%%% %%%%%%%%%%%%% Plot illustrative input-output curves %%%%%%%%%%%%%%%%%%%%

% % identify  systems near optimal sharpness for each gene circuit that
% % attain maximum near HM
% rate_bounds = [0.49 0.51];
% cw_val_vec = [1e1,1e2,1e3];
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % find optimal neq systems
% pd_vec_neq = metric_array_neq(:,rate_index);
% cw_vec_neq_raw = metric_array_neq(:,cw_index);
% cw_vec_neq = 10.^metric_array_neq(:,cw_index);
% rate_array_neq = vertcat(sim_struct_neq.rate_array);
% opt_indices = [];
% opt_vals = [];
% 
% pd_filter_neq = pd_vec_neq>=rate_bounds(1)&pd_vec_neq<=rate_bounds(2)&sharpness_vec_neq>=0;
% 
% % define vector of cr values
% cr_vec = logspace(-2,2,1e3)';
% ir_rate_array = NaN(length(cr_vec),size(rate_array_neq,2),length(cw_val_vec));
% 
% for c = 1:length(cw_val_vec)
%     cw_filter = cw_vec_neq >=cw_val_vec(c)*0.95 & cw_vec_neq <=cw_val_vec(c)*1.05;
%     temp_options = find(cw_filter&pd_filter_neq);
%     [ir_max,ir_i_neq] = nanmax(ir_vec_neq(temp_options));
%     opt_indices(end+1) = temp_options(ir_i_neq);
%     opt_vals(end+1) = ir_max;
%     
%     opt_rates = rate_array_neq(opt_indices(end),:);
%     temp_rates = repmat(opt_rates,length(cr_vec),1);
%     temp_rates(:,sim_info_neq.cr_index) = cr_vec';
%     ir_rate_array(:,:,c) = temp_rates;
% end    
% 
% % map to correct subfolder
% functionPath = sim_info_neq.functionPath;
% slashes = unique([strfind(functionPath,'\') strfind(functionPath,'/')]);
% simName = functionPath(slashes(end-1)+1:slashes(end)-1);
% rmpath(genpath('../utilities/metricFunctions/'));
% addpath(genpath(['../utilities/metricFunctions/' simName]));
% 
% % calculate production rates and noise profiles
% pd_array = NaN(length(cr_vec),length(cw_val_vec));
% noise_array = NaN(length(cr_vec),length(cw_val_vec));
% for c = 1:length(cw_val_vec)
%   
%     temp_rates = ir_rate_array(:,:,c);         
%          
% %     paramCellInitNeq = mat2cell(temp_rates,size(temp_rates,1),ones(1,size(temp_rates,2)));
% %     
% %     TauOn = TauONFunction(paramCellInitNeq{:});
% %     TauOff = TauOFFFunction(paramCellInitNeq{:});
% %     TauCycle = TauOff+TauOn;
% %     temp_rates(:,4:end) = temp_rates(:,4:end).*TauCycle;  
%     paramCellNeq = mat2cell(temp_rates,size(temp_rates,1),ones(1,size(temp_rates,2)));
%     % neq
%     pd_array(:,c) = productionRateFunction(paramCellNeq{:});
%     noise_array(:,c) = sqrt(intrinsicVarianceFunction(paramCellNeq{:}));    
%     
% %     TauOn = TauONFunction(paramCellNeq{:});
% %     TauOff = TauOFFFunction(paramCellNeq{:});
% %     TauCycle = TauOff+TauOn;
% end
% 
% 
% close all
% c1 = sim_info_neq.cr1;
% c0 = sim_info_neq.cr0;
% 
% 
% induction_plots_sym = figure;
% cmap = brewermap(8,'Blues');
% 
% hold on
% n_cycles = 50;%60/5;
% pneq = [];
% for c = length(cw_val_vec):-1:1
%   
%     % calculate upper and lower bounds 
%     ubneq = pd_array(:,c) + noise_array(:,c)/sqrt(n_cycles);
%     lbneq = pd_array(:,c) - noise_array(:,c)/sqrt(n_cycles);
% 
%     % plot verticle lines denoting target concentrations
%     plot(repelem(c0,100), linspace(0,1,100),'.k','LineWidth',1.5)
%     plot(repelem(c1,100), linspace(0,1,100),'.k','LineWidth',1.5)
% 
%     % make area plots to indicate error range for each    
%     fill([cr_vec' fliplr(cr_vec')],[lbneq' fliplr(ubneq')],cmap(c*2,:),'FaceAlpha',.2,'EdgeAlpha',0.75,'EdgeColor',cmap(2*c,:));
%     pneq(end+1) = plot(cr_vec,pd_array(:,c) ,'Color',brighten(cmap(c*2,:),-0.25),'LineWidth',3);
%         
% end
% 
% ylim([0 1])
% xlim([0.1 10])
% set(gca,'xscale','log')
% 
% legend(fliplr(pneq),'w/c = 10','w/c = 10^2','w/c = 10^3','Location','southeast')
% xlabel('activator concentration (c)');
% ylabel('transcription rate (r)')
% 
% grid on
% 
% set(gca,'FontSize',14)
% set(gca,'Color',[228,221,209]/255) 
% 
% induction_plots_sym.InvertHardcopy = 'off';
% set(gcf,'color','w');
% 
% saveas(induction_plots_sym,[FigPath 'cw_induction.png'])
% saveas(induction_plots_sym,[FigPath 'cw_induction.pdf'])    
% 
% 
% %% %%%%%%%%%%%%%%%% now look at decision time vs cw %%%%%%%%%%%%%%%%%%%%%%%
% close all
% 
% rate_array_neq = vertcat(sim_struct_neq.rate_array);
% rate_array_eq = vertcat(sim_struct_eq.rate_array);
% 
% % close all
% err = 0.32;
% t_cycle = 1;%5*60; % set the cycle time
% cw_vec = log10(logspace(0,log10(alpha_factor^4.1),1e2+1));
% 
% % initialize arrays to store rates
% tau_array = NaN(length(cw_vec)-1,2);
% 
% % iterate through cw values and find optimal eq and non-eq networks
% for c = 1:length(cw_vec)-1
%   % neq
%   cw_filter_neq = cw_vec_neq_raw>=cw_vec(c) & cw_vec_neq_raw<cw_vec(c+1);
%   [neq_max, optimal_neq_index] = nanmax(ir_vec_neq.*cw_filter_neq);
%   
%   % eq
%   cw_filter_eq = cw_vec_eq>=cw_vec(c) & cw_vec_eq<cw_vec(c+1);
%   [eq_max, optimal_eq_index] = nanmax(ir_vec_eq.*cw_filter_eq);
%   
%   % extract cycle times
%   % neq
%   cycle_time_neq = metric_array_neq(optimal_neq_index,cycle_time_index);
%   optimal_rates_neq = rate_array_neq(optimal_neq_index,:);%*cycle_time_neq/t_cycle;
%   optimal_rates_neq(4:end) = optimal_rates_neq(4:end)*cycle_time_neq/t_cycle;
%   
%   % eq
%   cycle_time_eq = metric_array_eq(optimal_eq_index,cycle_time_index);
%   optimal_rates_eq = rate_array_eq(optimal_eq_index,:);% % adjust to proper cycle time
%   optimal_rates_eq(4:end) = optimal_rates_eq(4:end)*cycle_time_eq/t_cycle;
% 
%   % calculate decision times
%   % neq
%   sim_info_neq_temp = sim_info_neq;
%   sim_info_neq_temp.minError = err;
%   [Veq,tau_array(c,2),~] = calculateDecisionMetrics(optimal_rates_neq,sim_info_neq_temp);     
%   
%   % eq
%   sim_info_eq_temp = sim_info_eq;
%   sim_info_eq_temp.minError = err;
%   [Vneq,tau_array(c,1),~] = calculateDecisionMetrics(optimal_rates_eq,sim_info_eq_temp);     
% end
% 
% 
% cw_vec_plot = 10.^(cw_vec(1:end-1) + diff(cw_vec)/2);
% cw_vec_plot(1) = 1;
% 
% [~,mi] = min(abs(cw_vec_plot-50))
% tau_array(mi,:)/t_cycle
% tau_array(mi,:)/60/60
% 
% 
% % generate decision time range vectors
% dt_cw_vec = geomean(vertcat(decision_limit_info.cw_ub,decision_limit_info.cw_lb));
% dt_ub_vec = decision_limit_info.n_cycles_ub;
% dt_lb_vec = max(vertcat(decision_limit_info.n_cycles_lb,ones(size(dt_cw_vec)))); % cap lower limit at 1 cycle
% dt_mean_vec = mean(vertcat(dt_ub_vec,dt_lb_vec));
% 
% % combine mous and human
% % dt_cw_vec(end-1) = mean(dt_cw_vec(end-1:end));
% % dt_cw_vec = dt_cw_vec(1:end-1);
% % dt_ub_vec(end-1) = max(dt_ub_vec(end-1:end));
% % dt_ub_vec = dt_ub_vec(1:end-1);
% % dt_lb_vec(end-1) = min(dt_lb_vec(end-1:end));
% % dt_lb_vec = dt_lb_vec(1:end-1);
% 
% % dt_mean_vec(end-1) = mean(dt_mean_vec(end-1:end));
% % dt_mean_vec = dt_mean_vec(1:end-1);
% 
% % exclude Arabadopsis for now
% plot_vec = [1 1 0 0 1]==1;
% 
% 
% close all
% 
% %%%%%%%%%%%%%%%
% % make figure
% decision_time_plot = figure;%('Position',[100 100 560*1.15  420]);
% cmap2 = brewermap([],'Set2');
% 
% topline = repelem(1e17,length(cw_vec)-1);
% hold on
% 
% % plot eq and noneq regimes
% fill([cw_vec_plot fliplr(cw_vec_plot)], [tau_array(:,1)'/t_cycle topline],cmap2(3,:),'FaceAlpha',0.5,'EdgeAlpha',0)
% fill([cw_vec_plot fliplr(cw_vec_plot)], [tau_array(:,2)'/t_cycle fliplr(tau_array(:,1)')/t_cycle],cmap2(2,:),'FaceAlpha',0.5,'EdgeAlpha',0)
% 
% p1 = plot(cw_vec_plot,tau_array(:,1)/t_cycle,'-','Color',brighten(cmap2(3,:),-.2),'LineWidth',3); 
% p2 = plot(cw_vec_plot,imgaussfilt(tau_array(:,2)/t_cycle,1),'Color',brighten(cmap2(2,:),-.2),'LineWidth',3); 
% % plot(error_vec,tau_array(:,n_lines)/60,'Color','k','LineWidth',2); 
% 
% % plot decision time ranges
% errorbar(dt_cw_vec(plot_vec),dt_mean_vec(plot_vec),dt_mean_vec(plot_vec)-dt_lb_vec(plot_vec),dt_ub_vec(plot_vec)-dt_mean_vec(plot_vec),'o','LineWidth',1.5','Color','k')
% 
% % legend([p1 p2],'equilibrium limit','nonequilibrium limit','Location','southeast')    
% xlabel('non-cognate factor concentration (w/c)');
% ylabel('decision time (burst cycles)');
% 
% grid on
% set(gca,'FontSize',14)
% set(gca,'Color',[228,221,209]/255) 
% 
% ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.XAxis(1).Color = 'k';
% 
% xlim([1 1e5])
% ylim([1 1e5])
% set(gca,'yscale','log')
% set(gca,'xscale','log')
% % set(gca,'ytick',[1 1e5 1e10])
% set(gca,'xtick',[1 alpha_factor alpha_factor^2 alpha_factor^3 alpha_factor^4]);%,'xticklabels',{'\alpha^0','\alpha^1','\alpha^{2}','\alpha^3','\alpha^4'})
% decision_time_plot.InvertHardcopy = 'off';
% set(gcf,'color','w');
% 
% saveas(decision_time_plot,[FigPath 'decision_time_vs_cw.png'])
% saveas(decision_time_plot,[FigPath 'decision_time_vs_cw.pdf'])

% %% calculate relative non-eq gain as a function of CW
% close all
% gain_fig = figure;
% peq = plot(10.^cw_axis,imgaussfilt(neq_bound_vec./eq_bound_vec,2),'Color','k','LineWidth',3);
% set(gca,'xscale','log') 
% set(gca,'yscale','log') 
% xlabel('non-cognate factor concentration (w/c)');
% ylabel('non-equilibrium gain (IR_{neq}/IR_{eq})')
% grid on
% set(gca,'FontSize',14)
% set(gca,'Color',[228,221,209]/255) 
% set(gca,'xtick',[1 alpha_factor alpha_factor^2 alpha_factor^3 alpha_factor^4],'xticklabels',{'\alpha^0','\alpha^1','\alpha^{2}','\alpha^3','\alpha^4'})
% set(gca,'ytick',[1  alpha_factor^1  alpha_factor^2],'yticklabels',{'\alpha^0','\alpha^{1}','\alpha^2'})
% ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.XAxis(1).Color = 'k';
% % ylim([10^-6.7 1])
% xlim([1 alpha_factor^4])
% grid on
% gain_fig.InvertHardcopy = 'off';
% set(gcf,'color','w');
% 
% saveas(gain_fig,[FigPath 'neq_ir_gain.png'])
% saveas(gain_fig,[FigPath 'neq_ir_gain.pdf'])