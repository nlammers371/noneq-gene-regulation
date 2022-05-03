% script to track information rate as a function of cw

clear 
close all
addpath(genpath('../utilities/'))

% set basic parameters
nStates = 6;
paramBounds = repmat([-10 ; 6],1,11); % constrain transition rate magnitude
% paramBounds(2,[8,10]) = 2;
[~,~,metric_names] = calculateMetricsSym_v2([]);

% specify function path
functionPath = ['../utilities/metricFunctions/n' num2str(nStates) '_OR/'];

% make sure we're linked to the appropriate function subfolder% make sure we're linked to the appropriate function subfolder
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));

% define save path

DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\manuscript\';
% DropboxFolder = 'S:\Nick\Dropbox\Nonequilibrium\Nick\manuscript\';

FigPath = [DropboxFolder 'information_vs_cw' filesep];
mkdir(FigPath);         

% get index of useful metrics
ir_index = find(strcmp(metric_names,'DecisionRateNorm'));
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
TauCycleLimit = 300; % do this for circuits with a 5 minute cycle time

bit_factor = log2(exp(1));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% info vs cw %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [sim_info_neq, sim_struct_neq] = param_sweep_multi_v3([ir_index phi_index],functionPath,sweep_options{:},...
%                                           'half_max_flag',false,'TauCycleLimit',TauCycleLimit,...
%                                           'equilibrium_flag',false,'specFactor',alpha_factor,'paramBounds',paramBounds,'cw',1e-6);
% 
% %%
tic
[sim_info_neq, sim_struct_neq] = param_sweep_multi_v3([ir_index cw_index],functionPath,sweep_options{:},...
                                          'half_max_flag',false,'TauCycleLimit',TauCycleLimit,...
                                          'equilibrium_flag',false,'specFactor',alpha_factor,'paramBounds',paramBounds);

[sim_info_eq, sim_struct_eq] = param_sweep_multi_v3([ir_index cw_index],functionPath,sweep_options{:},...
                                          'half_max_flag',false,'TauCycleLimit',TauCycleLimit,...
                                          'equilibrium_flag',true,'specFactor',alpha_factor,'paramBounds',paramBounds);                                        
toc     

%% pull result arrays
n_plot = 3e3;
metric_array_neq = vertcat(sim_struct_neq.metric_array);
% rate_array_neq = vertcat(sim_struct_neq.rate_array);
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
cw_vec = log10(logspace(log10(alpha_factor/1000),log10(alpha_factor^4.1),31));
neq_bound_vec = NaN(1,length(cw_vec)-1);
eq_bound_vec = NaN(1,length(cw_vec)-1);

for c = 1:length(cw_vec)-1
    cw_neq_filter = cw_vec_neq(neq_indices)>=cw_vec(c) & cw_vec_neq(neq_indices)<cw_vec(c+1);
    if any(cw_neq_filter)
      neq_bound_vec(c) = nanmax(ir_vec_neq(neq_indices(cw_neq_filter)));
    end
    cw_eq_filter = cw_vec_eq(eq_indices)>=cw_vec(c) & cw_vec_eq(eq_indices)<cw_vec(c+1);
    if any(cw_eq_filter)
      eq_bound_vec(c) = nanmax(ir_vec_eq(eq_indices(cw_eq_filter)));
    end
end
cw_axis = cw_vec(1:end-1);
% cw_axis(1) = 0;

%% make figure
ir_vs_cw_eq_fig = figure;
hold on
cmap = brewermap([],'Set2');

sneq = scatter(10.^cw_vec_neq(plot_indices_neq), ir_vec_neq(plot_indices_neq),...
      markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',0.25, 'MarkerFaceColor',cmap(2,:));
    
seq = scatter(10.^cw_vec_eq(plot_indices_eq), ir_vec_eq(plot_indices_eq),...
      markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',0.25, 'MarkerFaceColor',cmap(3,:));    
    
plot(10.^cw_axis,imgaussfilt(neq_bound_vec),'--','Color',brighten(cmap(2,:),-.5),'LineWidth',3)
plot(10.^cw_axis,imgaussfilt(eq_bound_vec),'--','Color',brighten(cmap(3,:),-.5),'LineWidth',3)

set(gca,'xscale','log') 
set(gca,'yscale','log') 
xlabel('non-cognate factor concentration (c_w/c_r)');
ylabel('information rate (bits per cycle)')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
set(gca,'xtick',[1 alpha_factor alpha_factor^2 alpha_factor^3 alpha_factor^4],'xticklabels',{'\alpha^0','\alpha^1','\alpha^{2}','\alpha^3','\alpha^4'})
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
ylim([1e-10 1])
xlim([alpha_factor/1e3 alpha_factor^4])
grid on
ir_vs_cw_eq_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
% saveas(ir_vs_cw_eq_fig,[FigPath 'IR_vs_cw.png'])
% saveas(ir_vs_cw_eq_fig,[FigPath 'IR_vs_cw.pdf'])



%% calculate relative non-eq gain as a function of CW


close all
gain_fig = figure;
peq = plot(10.^cw_axis,imgaussfilt(neq_bound_vec./eq_bound_vec,2),'Color','k','LineWidth',3);
set(gca,'xscale','log') 
set(gca,'yscale','log') 
xlabel('non-cognate factor concentration (c_w/c_r)');
ylabel('non-equilibrium gain (IR_{neq}/IR_{eq})')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
set(gca,'xtick',[1 alpha_factor alpha_factor^2 alpha_factor^3 alpha_factor^4],'xticklabels',{'\alpha^0','\alpha^1','\alpha^{2}','\alpha^3','\alpha^4'})
set(gca,'ytick',[1  alpha_factor^1  alpha_factor^2],'yticklabels',{'\alpha^0','\alpha^{1}','\alpha^2'})
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
% ylim([10^-6.7 1])
xlim([1 alpha_factor^4])
grid on
gain_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(gain_fig,[FigPath 'neq_ir_gain.png'])
saveas(gain_fig,[FigPath 'neq_ir_gain.pdf'])

%% %%%%%%%%%%%%%%%% now look at decision time vs cw %%%%%%%%%%%%%%%%%%%%%
close all

rate_array_neq = vertcat(sim_struct_neq.rate_array);
rate_array_eq = vertcat(sim_struct_eq.rate_array);

% close all
err = 0.32;
t_cycle = 300;%5*60; % set the cycle time
cw_vec = log10(logspace(0,log10(alpha_factor^4.1),1e2+1));

% initialize arrays to store rates
tau_array = NaN(length(cw_vec)-1,2);

% iterate through cw values and find optimal eq and non-eq networks
for c = 1:length(cw_vec)-1
  % neq
  cw_filter_neq = cw_vec_neq>=cw_vec(c) & cw_vec_neq<cw_vec(c+1);
  [neq_max, optimal_neq_index] = nanmax(ir_vec_neq.*cw_filter_neq);
  
  % eq
  cw_filter_eq = cw_vec_eq>=cw_vec(c) & cw_vec_eq<cw_vec(c+1);
  [eq_max, optimal_eq_index] = nanmax(ir_vec_eq.*cw_filter_eq);
  
  % extract cycle times
  % neq
  cycle_time_neq = metric_array_neq(optimal_neq_index,cycle_time_index);
  optimal_rates_neq = rate_array_neq(optimal_neq_index,:);%*cycle_time_neq/t_cycle;
  optimal_rates_neq(4:end) = optimal_rates_neq(4:end)*cycle_time_neq/t_cycle;
  
  % eq
  cycle_time_eq = metric_array_eq(optimal_eq_index,cycle_time_index);
  optimal_rates_eq = rate_array_eq(optimal_eq_index,:);% % adjust to proper cycle time
  optimal_rates_eq(4:end) = optimal_rates_eq(4:end)*cycle_time_eq/t_cycle;

  % calculate decision times
  % neq
  sim_info_neq_temp = sim_info_neq;
  sim_info_neq_temp.minError = err;
  [Veq,tau_array(c,2),~] = calculateDecisionMetrics(optimal_rates_neq,sim_info_neq_temp);     
  
  % eq
  sim_info_eq_temp = sim_info_eq;
  sim_info_eq_temp.minError = err;
  [Vneq,tau_array(c,1),~] = calculateDecisionMetrics(optimal_rates_eq,sim_info_eq_temp);     
end


cw_vec_plot = 10.^(cw_vec(1:end-1) + diff(cw_vec)/2);
cw_vec_plot(1) = 1;

[~,mi] = min(abs(cw_vec_plot-alpha_factor^2))
tau_array(mi,:)/t_cycle
tau_array(mi,:)/60/60

%%
close all
decision_time_plot = figure;%('Position',[100 100 560*1.15  420]);
% cmap = flipud(brewermap(n_lines,'Spectral'));
cmap2 = brewermap([],'Set2');


% cw_vec_plot(end) = 1e6;

topline = repelem(1e17,length(cw_vec)-1);
hold on

% plot eq and noneq regimes
fill([cw_vec_plot fliplr(cw_vec_plot)], [tau_array(:,1)'/t_cycle topline],cmap2(3,:),'FaceAlpha',0.5,'EdgeAlpha',0)
fill([cw_vec_plot fliplr(cw_vec_plot)], [tau_array(:,2)'/t_cycle fliplr(tau_array(:,1)')/t_cycle],cmap2(2,:),'FaceAlpha',0.5,'EdgeAlpha',0)

p1 = plot(cw_vec_plot,tau_array(:,1)/t_cycle,'-','Color',brighten(cmap2(3,:),-.2),'LineWidth',3); 
p2 = plot(cw_vec_plot,imgaussfilt(tau_array(:,2)/t_cycle,1),'Color',brighten(cmap2(2,:),-.2),'LineWidth',3); 
% plot(error_vec,tau_array(:,n_lines)/60,'Color','k','LineWidth',2); 

legend([p1 p2],'equilibrium limit','nonequilibrium limit','Location','southeast')    
xlabel('non-cognate factor concentration (c_w/c_r)');
ylabel('\langle T \rangle (burst cycles)');


grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

xlim([1 alpha_factor^4])
ylim([1 1e10])
set(gca,'yscale','log')
set(gca,'xscale','log')
set(gca,'ytick',[1 1e5 1e10])
set(gca,'xtick',[1 alpha_factor alpha_factor^2 alpha_factor^3 alpha_factor^4],'xticklabels',{'\alpha^0','\alpha^1','\alpha^{2}','\alpha^3','\alpha^4'})
decision_time_plot.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(decision_time_plot,[FigPath 'decision_time_vs_cw.png'])
saveas(decision_time_plot,[FigPath 'decision_time_vs_cw.pdf'])

% %% %%%%%%%%%%%%%%%% make phase space plots %%%%%%%%%%%%%%%%%%%%%
% close all
% 
% cw_vec1 = log10(logspace(0,6,50+1));
% 
% precision_vec_neq = exp(metric_array_neq(:,precision_index)).^2;
% spec_vec_neq = 10.^metric_array_neq(:,spec_index);
% s0_vec_neq = metric_array_neq(:,sharp_right_norm_index);
% 
% p_max_eq = exp(nanmax(metric_array_eq(:,precision_index))).^2;
% f_max_eq = 10.^nanmax(metric_array_eq(:,spec_index));
% s_max_eq = nanmax(metric_array_eq(:,sharp_right_norm_index));
% 
% ir_vec = NaN(length(cw_vec1)-1,1);
% p_vec = NaN(length(cw_vec1)-1,1);
% s_vec = NaN(length(cw_vec1)-1,1);
% f_vec = NaN(length(cw_vec1)-1,1);
% 
% for c = 1:length(cw_vec1)-1
%   % neq
%   cw_filter_neq = cw_vec_neq>=cw_vec1(c) & cw_vec_neq<cw_vec1(c+1) & s0_vec_neq>=0.99;
%   [ir_vec(c), optimal_neq_index] = nanmax(ir_vec_neq.*cw_filter_neq);    
%   
%   % store metric values for optimal network
%   % neq
%   p_vec(c) = precision_vec_neq(optimal_neq_index)/p_max_eq;
%   s_vec(c) = s0_vec_neq(optimal_neq_index)/s_max_eq;
%   f_vec(c) = spec_vec_neq(optimal_neq_index)/f_max_eq;
% end
% 
% cw_vec1_plot = 10.^(cw_vec1(1:end-1) + diff(cw_vec1)/2);
% 
% % sharpness ves specificity
% sharp_vs_spec = figure;
% cmap2 = flipud(brewermap(length(cw_vec1)-1,'Spectral'));
% colormap(cmap2)
% hold on
% 
% plot(linspace(0.01,alpha_factor.^2),ones(1,100),'--k','LineWidth',2)
% plot(ones(1,100),linspace(0,4),'--k','LineWidth',2)
% 
% for c = 1:length(cw_vec1_plot)-1
%     seq = scatter(f_vec(c).^2, s_vec(c).^2,...
%               markerSize*1.5,'MarkerEdgeAlpha',.5,'MarkerEdgeColor','k',...
%               'MarkerFaceAlpha',0.75, 'MarkerFaceColor',cmap2(c,:));      
% end  
% 
% xlim([.01 alpha_factor^2])
% ylim([0 4])
% 
% ylabel('sharpness gain (s^2/s^2_{eq})');
% % xlabel('precision gain (\sigma^2_{eq}/\sigma^2)')
% xlabel('specificity gain (f^2/f_{eq}^2)')
% 
% grid on
% set(gca,'FontSize',14)
% set(gca,'Color',[228,221,209]/255) 
% 
% ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.XAxis(1).Color = 'k';
% 
% h = colorbar;
% ylabel(h,'c_w/c_r')
% set(h,'xTickLabel',{'\alpha^0','\alpha^1','\alpha^2','\alpha^3'},...
%                'xTick', linspace(0,1,4))
% 
% % set(gca,'yscale','log')
% set(gca,'xscale','log')
% set(gca,'xtick',[alpha_factor.^-1 1 alpha_factor.^1 alpha_factor^2],'xticklabels',{'\alpha^{-1}','\alpha^0','\alpha^{1}','\alpha^{2}'})
% sharp_vs_spec.InvertHardcopy = 'off';
% set(gcf,'color','w');
% 
% saveas(sharp_vs_spec,[FigPath 'sharp_vs_spec.png'])
% saveas(sharp_vs_spec,[FigPath 'sharp_vs_spec.pdf'])
% 
% 
% % precision vs specificity
% 
% prec_vs_spec = figure;
% colormap(cmap2)
% hold on
% 
% plot(linspace(0.01,alpha_factor.^2),ones(1,100),'--k','LineWidth',2)
% plot(ones(1,100),linspace(0,16),'--k','LineWidth',2)
% 
% for c = 1:length(cw_vec1_plot)-1
%     seq = scatter(f_vec(c).^2, p_vec(c),...
%               markerSize*1.5,'MarkerEdgeAlpha',.5,'MarkerEdgeColor','k',...
%               'MarkerFaceAlpha',0.75, 'MarkerFaceColor',cmap2(c,:));      
% end  
% 
% xlim([.01 alpha_factor^2])
% ylim([0 2])
% 
% % ylabel('sharpness gain (s^2/s^2_{eq})');
% ylabel('precision gain (\sigma^2_{eq}/\sigma^2)')
% xlabel('specificity gain (f^2/f_{eq}^2)')
% 
% grid on
% set(gca,'FontSize',14)
% set(gca,'Color',[228,221,209]/255) 
% 
% ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.XAxis(1).Color = 'k';
% 
% h = colorbar;
% ylabel(h,'c_w/c_r')
% set(h,'xTickLabel',{'\alpha^0','\alpha^1','\alpha^2','\alpha^3'},...
%                'xTick', linspace(0,1,4))
% 
% % set(gca,'yscale','log')
% set(gca,'xscale','log')
% set(gca,'xtick',[alpha_factor.^-1 1 alpha_factor.^1 alpha_factor^2],'xticklabels',{'\alpha^{-1}','\alpha^0','\alpha^{1}','\alpha^{2}'})
% prec_vs_spec.InvertHardcopy = 'off';
% set(gcf,'color','w');
% 
% saveas(prec_vs_spec,[FigPath 'prec_vs_spec.png'])
% saveas(prec_vs_spec,[FigPath 'prec_vs_spec.pdf'])
% 
% 
% prec_vs_sharp = figure;
% colormap(cmap2)
% hold on
% 
% plot(linspace(0,2),ones(1,100),'--k','LineWidth',2)
% plot(ones(1,100),linspace(0,4),'--k','LineWidth',2)
% 
% for c = 1:length(cw_vec1_plot)-1
%     seq = scatter(p_vec(c), s_vec(c).^2,...
%               markerSize*1.5,'MarkerEdgeAlpha',.5,'MarkerEdgeColor','k',...
%               'MarkerFaceAlpha',0.75, 'MarkerFaceColor',cmap2(c,:));      
% end  
% 
% xlim([0 2])
% ylim([0 4])
% 
% ylabel('sharpness gain (s^2/s^2_{eq})');
% xlabel('precision gain (\sigma^2_{eq}/\sigma^2)')
% % xlabel('specificity gain (f^2/f_{eq}^2)')
% 
% grid on
% set(gca,'FontSize',14)
% set(gca,'Color',[228,221,209]/255) 
% 
% ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.XAxis(1).Color = 'k';
% 
% h = colorbar;
% ylabel(h,'c_w/c_r')
% set(h,'xTickLabel',{'\alpha^0','\alpha^1','\alpha^2','\alpha^3'},...
%                'xTick', linspace(0,1,4))
% 
% % set(gca,'yscale','log')
% % set(gca,'xscale','log')
% % set(gca,'xtick',[alpha_factor.^-1 1 alpha_factor.^1 alpha_factor^2],'xticklabels',{'\alpha^{-1}','\alpha^0','\alpha^{1}','\alpha^{2}'})
% prec_vs_sharp.InvertHardcopy = 'off';
% set(gcf,'color','w');
% 
% saveas(prec_vs_sharp,[FigPath 'prec_vs_sharp.png'])
% saveas(prec_vs_sharp,[FigPath 'prec_vs_sharp.pdf'])