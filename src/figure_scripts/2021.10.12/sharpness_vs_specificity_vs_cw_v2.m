% script to track information rate as a function of cw

clear 
close all
addpath(genpath('../utilities/'))

% set basic parameters
nStates = 6;
paramBounds = repmat([-10 ; 6],1,3*nStates-4); % constrain transition rate magnitude
[~,~,metric_names] = calculateMetricsSym_v2([]);

% specify function path
functionPath = ['../utilities/metricFunctions/n' num2str(nStates) '_OR/'];

% make sure we're linked to the appropriate function subfolder% make sure we're linked to the appropriate function subfolder
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));

% define save path

DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\manuscript\';
% DropboxFolder = 'S:\Nick\Dropbox\Nonequilibrium\Nick\manuscript\';

FigPath = [DropboxFolder 'f0_s0_ir_cw' filesep];
mkdir(FigPath);         

% get index of useful metrics
spec_index = find(strcmp(metric_names,'Specificity'));
sharp_right_norm_index = find(strcmp(metric_names,'SharpnessRightNorm'));
ir_index = find(strcmp(metric_names,'DecisionRateNorm'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
precision_index = find(strcmp(metric_names,'Precision'));
precision_right_index = find(strcmp(metric_names,'PrecisionRight'));
cw_index = find(strcmp(metric_names,'CW'));

% set sim options
sweep_options = {'n_seeds',5,'n_iters_max',50,'n_sim',10,'nStates',nStates};

% calculate sensitivity bound
alpha_factor = 100;
% seq = 1/4;

% specify plot options 
markerAlpha = 0.5; % marker transparency
markerSize = 75; % marker size
n_plot = 3e3;
bit_factor = log2(exp(1));

%% Plot f0 vs s0 
% cw = 1;
% [sim_info_neq, sim_struct_neq] = param_sweep_multi_v3([precision_index spec_index],functionPath,sweep_options{:},...
%                                           'half_max_flag',true,'cw',cw,...
%                                           'equilibrium_flag',false,'specFactor',alpha_factor);
% 
% %%
cw = 1; % precise value unimportant...just something small enough to be negligible
tic
[sim_info_neq, sim_struct_neq] = param_sweep_multi_v3([sharp_right_norm_index spec_index],functionPath,sweep_options{:},...
                                          'half_max_flag',false,'cw',cw,...
                                          'equilibrium_flag',false,'specFactor',alpha_factor);

[sim_info_eq, sim_struct_eq] = param_sweep_multi_v3([sharp_right_norm_index spec_index],functionPath,sweep_options{:},...
                                          'half_max_flag',false,'cw',cw,...
                                          'equilibrium_flag',true,'specFactor',alpha_factor);                                        
toc     

close all

% cw = 1;
% get predicted tradeoff bound
% s0_bound = seq+(alpha_factor.^2+(-1).*f0_vec).*cw.*((-1).*alpha_factor+f0_vec+((-1)+alpha_factor).*f0_vec.*cw).^( ...
%   -1).*seq;
seq = 1/4;
cr = 1;
n_plot1 = 3e3;
% s0_bound = (cr.^(-1)+((-1)+alpha_factor).^(-1).*cr.^(-1).*(alpha_factor.^2+(-1).*f0_vec).*f0_vec.^(-1));             
% s0_bound2 = (cr.^(-1)+(alpha_factor.^2+(-1).*f0_vec).*((-1).*alpha_factor.*cw+(((-1)+alpha_factor).*cr+cw).*f0_vec).^(-1));

% generate vectors to plot (neq)
results_array_neq = vertcat(sim_struct_neq.metric_array);
f0_scatter_vec_neq = 10.^results_array_neq(:,spec_index);
s0_scatter_vec_neq = results_array_neq(:,sharp_right_norm_index);
plot_options_neq = find(s0_scatter_vec_neq>=0 & f0_scatter_vec_neq >= 1);
plot_indices_neq = randsample(plot_options_neq,min([n_plot1,length(plot_options_neq)]),false);

% generate vectors to plot (eq)
results_array_eq = vertcat(sim_struct_eq.metric_array);
f0_scatter_vec_eq = 10.^results_array_eq(:,spec_index);
s0_scatter_vec_eq = results_array_eq(:,sharp_right_norm_index);
plot_options_eq = find(s0_scatter_vec_eq>=0 & f0_scatter_vec_eq >= 1);
plot_indices_eq = randsample(plot_options_eq,min([n_plot1,length(plot_options_eq)]),false);

% make figure
s0_f0_fig = figure;
hold on
cmap = brewermap([],'Set2');

% sneq = scatter(f0_scatter_vec_neq(plot_indices_neq)*alpha_factor,s0_scatter_vec_neq(plot_indices_neq),...
%       markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha*0.5, 'MarkerFaceColor',cmap(2,:));
    
seq = scatter(f0_scatter_vec_eq(plot_indices_eq)*alpha_factor,s0_scatter_vec_eq(plot_indices_eq),...
      markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha*0.5, 'MarkerFaceColor',cmap(3,:));    
    
% p = plot(f0_vec,s0_bound2,'-.','Color','k','LineWidth',3);
ylim([0 2])
xlim([alpha_factor alpha_factor^2])
set(gca,'xscale','log')    
xlabel('intrinsic specificity (f_0)');
ylabel('intrinsic sharpness (s_0/s_0^{eq})')
% grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
set(gca,'xtick',[alpha_factor alpha_factor^1.5 alpha_factor^2],'xticklabels',{'\alpha','\alpha^{1.5}','\alpha^2'})
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

s0_f0_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
% set(p,'Visible','off')

saveas(s0_f0_fig,[FigPath 's0_vs_f0.png'])

% now add non-equilibrium results
sneq = scatter(f0_scatter_vec_neq(plot_indices_neq)*alpha_factor,s0_scatter_vec_neq(plot_indices_neq),...
      markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha*0.5, 'MarkerFaceColor',cmap(2,:));
uistack(seq,'top')
% saveas(s0_f0_fig,[FigPath 's0_vs_f0_neq.png'])
% saveas(s0_f0_fig,[FigPath 's0_vs_f0_neq.pdf'])


%% %%%%%%%%%%%%%%%%%%%%%%%%%%% info vs cw %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
[sim_info_neq, sim_struct_neq] = param_sweep_multi_v3([ir_index cw_index],functionPath,sweep_options{:},...
                                          'half_max_flag',true,...
                                          'equilibrium_flag',false,'specFactor',alpha_factor,'paramBounds',paramBounds);

[sim_info_eq, sim_struct_eq] = param_sweep_multi_v3([ir_index cw_index],functionPath,sweep_options{:},...
                                          'half_max_flag',true,...
                                          'equilibrium_flag',true,'specFactor',alpha_factor,'paramBounds',paramBounds);                                        
toc     

% pull result arrays
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

%%%%%%%%%%%%%%%%%%%%
% Now for s0 and f0-maximized motifs

tic
[sim_info_s0, sim_struct_s0] = param_sweep_multi_v3([sharp_right_norm_index cw_index],functionPath,sweep_options{:},...
                                          'half_max_flag',true,...
                                          'equilibrium_flag',false,'specFactor',alpha_factor,'paramBounds',paramBounds);

[sim_info_f0, sim_struct_f0] = param_sweep_multi_v3([spec_index cw_index],functionPath,sweep_options{:},...
                                          'half_max_flag',true,...
                                          'equilibrium_flag',false,'specFactor',alpha_factor,'paramBounds',paramBounds);                                        
toc   
close all


%% pull result arrays
metric_array_f0 = vertcat(sim_struct_f0.metric_array);
metric_array_s0 = vertcat(sim_struct_s0.metric_array);

% extract key vectors
s0_vec_f0 = metric_array_f0(:,sharp_right_norm_index);
s0_vec_s0 = metric_array_s0(:,sharp_right_norm_index);
f0_vec_s0 = metric_array_s0(:,spec_index);
f0_vec_f0 = metric_array_f0(:,spec_index);

sharpness_vec_f0 = metric_array_f0(:,sharpness_index);
sharpness_vec_s0 = metric_array_s0(:,sharpness_index);

ir_vec_f0 = metric_array_f0(:,ir_index)*bit_factor;
ir_vec_s0 = metric_array_s0(:,ir_index)*bit_factor;

cw_vec_f0 = metric_array_f0(:,cw_index);
cw_vec_s0 = metric_array_s0(:,cw_index);

% create filters
s0_indices = find(f0_vec_s0<=1 & s0_vec_s0 > 1.9);
f0_indices = find(s0_vec_f0<=1 & f0_vec_f0 > 1.9);

%% find upper IR limits for global, s0 and f0 regions
cw_bins = 10.^linspace(0,log10(alpha_factor^4),41);

ir_bound_f0 = NaN(1,length(cw_bins)-1);
ir_bound_s0 = NaN(1,length(cw_bins)-1);
ir_bound_global = NaN(1,length(cw_bins)-1);
ir_bound_global_eq = NaN(1,length(cw_bins)-1);

for c = 1:length(cw_bins)-1
    % f0 
    cw_filter_f0 = 10.^cw_vec_f0(f0_indices) >= cw_bins(c) & 10.^cw_vec_f0(f0_indices) < cw_bins(c+1);
    ir_bound_f0(c) = nanmax(ir_vec_f0(f0_indices(cw_filter_f0)));
    % s0 
    cw_filter_s0 = 10.^cw_vec_s0(s0_indices) >= cw_bins(c) & 10.^cw_vec_s0(s0_indices) < cw_bins(c+1);
    ir_bound_s0(c) = nanmax(ir_vec_s0(s0_indices(cw_filter_s0)));
    % global 
    cw_filter_neq = 10.^cw_vec_neq(neq_indices) >= cw_bins(c) & 10.^cw_vec_neq(neq_indices) < cw_bins(c+1);
    ir_bound_global(c) = nanmax(ir_vec_neq(neq_indices(cw_filter_neq)));
    % global eq
    cw_filter_eq = 10.^cw_vec_eq(eq_indices) >= cw_bins(c) & 10.^cw_vec_eq(eq_indices) < cw_bins(c+1);
    ir_bound_global_eq(c) = nanmax(ir_vec_eq(eq_indices(cw_filter_eq)));
    
end
% add replicate first point 
ir_bound_f0 = [ir_bound_f0(1) ir_bound_f0];
ir_bound_s0 = [ir_bound_s0(1) ir_bound_s0];
ir_bound_global = [ir_bound_global(1) ir_bound_global];
ir_bound_global_eq = [ir_bound_global_eq(1) ir_bound_global_eq];

% Make plots
close all

% make figure
ir_vs_cw_eq_fig = figure;
hold on
cmap = brewermap([],'Set2');

pneq = plot(cw_bins, ir_bound_global, 'Color',cmap(2,:), 'LineWidth', 3);
    
peq = plot(cw_bins, ir_bound_global_eq, 'Color',cmap(3,:), 'LineWidth', 3);

ps0 = scatter(cw_bins, ir_bound_s0, 50,'s', 'MarkerFaceColor',cmap(4,:),'MarkerEdgeColor','k');

pf0 = scatter(cw_bins, ir_bound_f0, 50,'^', 'MarkerFaceColor',cmap(5,:),'MarkerEdgeColor','k');
    
% ylim([0 2])

set(gca,'xscale','log') 
set(gca,'yscale','log') 
xlabel('non-cognate factor concentration (c_w/c_r)');
ylabel('information rate (bits per cycle)')
grid on
set(gca,'FontSize',14)

legend([pneq peq ps0 pf0],'non-eq. optimum','eq. optimum','s_0-optimized','f_0-optimized','Location','southwest')

set(gca,'Color',[228,221,209]/255) 
set(gca,'xtick',[1 alpha_factor alpha_factor^2 alpha_factor^3 alpha_factor^4],'xticklabels',{'\alpha^0','\alpha^1','\alpha^{2}','\alpha^3','\alpha^4'})
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
ylim([1e-10 1])
xlim([1 alpha_factor^4])
grid on
ir_vs_cw_eq_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
% saveas(ir_vs_cw_eq_fig,[FigPath 'ir_s0_vs_f0_vs_global.png'])
% saveas(ir_vs_cw_eq_fig,[FigPath 'ir_s0_vs_f0_vs_global.pdf'])


%% %%%%%%%%%%%%%%%% make phase space plots %%%%%%%%%%%%%%%%%%%%%
close all
n_bins = 151;

cw_vec1 = log10(logspace(0,log10(alpha_factor^4),round(n_bins/8)*8+1));

precision_vec_neq = exp(metric_array_neq(:,precision_index)).^2;
spec_vec_neq = 10.^metric_array_neq(:,spec_index);
s0_vec_neq = metric_array_neq(:,sharp_right_norm_index);

p_max_eq = exp(nanmax(metric_array_eq(:,precision_index))).^2;
f_max_eq = 10.^nanmax(metric_array_eq(:,spec_index));
s_max_eq = nanmax(metric_array_eq(:,sharp_right_norm_index));

ir_vec = NaN(length(cw_vec1)-1,1);
p_vec = NaN(length(cw_vec1)-1,1);
s_vec = NaN(length(cw_vec1)-1,1);
f_vec = NaN(length(cw_vec1)-1,1);

for c = 1:length(cw_vec1)-1
  % neq
  cw_filter_neq = cw_vec_neq>=cw_vec1(c) & cw_vec_neq<cw_vec1(c+1) & s0_vec_neq>=0.99;
  [ir_vec(c), optimal_neq_index] = nanmax(ir_vec_neq.*cw_filter_neq);    
  
  % store metric values for optimal network
  % neq
  p_vec(c) = precision_vec_neq(optimal_neq_index)/p_max_eq;
  s_vec(c) = s0_vec_neq(optimal_neq_index)/s_max_eq;
  f_vec(c) = spec_vec_neq(optimal_neq_index)/f_max_eq;
end

cw_vec1_plot = 10.^(cw_vec1(1:end-1) + diff(cw_vec1)/2);

% generate raw colormap
c_unit = (length(cw_vec1)-1)*0.125;
cmap2 = flipud(brewermap(8*c_unit,'Spectral'));

% cmap_top1 = interp1(1:3*c_unit,cmap2(1:3*c_unit,:),linspace(1,3*c_unit,4*c_unit));
% cmap_bottom1 = interp1(1:5*c_unit,cmap2(3*c_unit+1:end,:),linspace(1,5*c_unit,4*c_unit));
% cmap3 = vertcat(cmap_top1,cmap_bottom1);
% 
cmap_bottom2 = interp1(1:3*c_unit,cmap2(5*c_unit+1:end,:),linspace(1,3*c_unit,4*c_unit));
cmap_top2 = interp1(1:5*c_unit,cmap2(1:5*c_unit,:),linspace(1,5*c_unit,4*c_unit));
cmap4 = vertcat(cmap_top2,cmap_bottom2);
% cmap4 = cmap2;

% sharpness ves specificity
sharp_vs_spec = figure;
colormap(cmap4)
hold on

plot(linspace(0.01,alpha_factor.^2),ones(1,100),'--k','LineWidth',2)
plot(ones(1,100),linspace(0,4),'--k','LineWidth',2)

for c = 1:length(cw_vec1_plot)-1
    peq = scatter(f_vec(c), s_vec(c),...
              markerSize*1.5,'MarkerEdgeAlpha',.5,'MarkerEdgeColor','k',...
              'MarkerFaceAlpha',0.75, 'MarkerFaceColor',cmap4(c,:));      
end  

xlim([.1 alpha_factor^1])
ylim([0.5 2])

ylabel('sharpness gain (s_0/s_0^{eq})');
% xlabel('precision gain (\sigma^2_{eq}/\sigma^2)')
xlabel('specificity gain (f_0/f_0^{eq})')

grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

h = colorbar;
ylabel(h,'c_w/c_r')
set(h,'xTickLabel',{'\alpha^0','\alpha^1','\alpha^2','\alpha^3','\alpha^4'},...
               'xTick', linspace(0,1,5))

% set(gca,'yscale','log')
set(gca,'xscale','log')
set(gca,'xtick',[alpha_factor.^-1 1 alpha_factor.^1 alpha_factor^2],'xticklabels',{'\alpha^{-1}','\alpha^0','\alpha^{1}','\alpha^{2}'})
sharp_vs_spec.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(sharp_vs_spec,[FigPath 'sharp_vs_spec.png'])
saveas(sharp_vs_spec,[FigPath 'sharp_vs_spec.pdf'])


%% precision vs specificity

prec_vs_spec = figure;
colormap(cmap4)
hold on

plot(linspace(0.01,alpha_factor.^2),ones(1,100),'--k','LineWidth',2)
plot(ones(1,100),linspace(0,16),'--k','LineWidth',2)

for c = 1:length(cw_vec1_plot)-1
    peq = scatter(f_vec(c).^2, p_vec(c),...
              markerSize*1.5,'MarkerEdgeAlpha',.5,'MarkerEdgeColor','k',...
              'MarkerFaceAlpha',0.75, 'MarkerFaceColor',cmap4(c,:));      
end  

xlim([.01 alpha_factor^2])
ylim([0 2])

% ylabel('sharpness gain (s^2/s^2_{eq})');
ylabel('precision gain (\sigma^2_{eq}/\sigma^2)')
xlabel('specificity gain (f^2/f_{eq}^2)')

grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

h = colorbar;
ylabel(h,'c_w/c_r')
set(h,'xTickLabel',{'\alpha^0','\alpha^1','\alpha^2','\alpha^3','\alpha^4'},...
               'xTick', linspace(0,1,5))

% set(gca,'yscale','log')
set(gca,'xscale','log')
set(gca,'xtick',[alpha_factor.^-1 1 alpha_factor.^1 alpha_factor^2],'xticklabels',{'\alpha^{-1}','\alpha^0','\alpha^{1}','\alpha^{2}'})
prec_vs_spec.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(prec_vs_spec,[FigPath 'prec_vs_spec.png'])
saveas(prec_vs_spec,[FigPath 'prec_vs_spec.pdf'])


prec_vs_sharp = figure;
colormap(cmap4)
hold on

plot(linspace(0,2),ones(1,100),'--k','LineWidth',2)
plot(ones(1,100),linspace(0,4),'--k','LineWidth',2)

for c = 1:length(cw_vec1_plot)-1
    peq = scatter(p_vec(c), s_vec(c).^2,...
              markerSize*1.5,'MarkerEdgeAlpha',.5,'MarkerEdgeColor','k',...
              'MarkerFaceAlpha',0.75, 'MarkerFaceColor',cmap4(c,:));      
end  

xlim([0 2])
ylim([0 4])

ylabel('sharpness gain (s^2/s^2_{eq})');
xlabel('precision gain (\sigma^2_{eq}/\sigma^2)')
% xlabel('specificity gain (f^2/f_{eq}^2)')

grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

h = colorbar;
ylabel(h,'c_w/c_r')
set(h,'xTickLabel',{'\alpha^0','\alpha^1','\alpha^2','\alpha^3','\alpha^4'},...
               'xTick', linspace(0,1,5))
% set(gca,'yscale','log')
% set(gca,'xscale','log')
% set(gca,'xtick',[alpha_factor.^-1 1 alpha_factor.^1 alpha_factor^2],'xticklabels',{'\alpha^{-1}','\alpha^0','\alpha^{1}','\alpha^{2}'})
prec_vs_sharp.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(prec_vs_sharp,[FigPath 'prec_vs_sharp.png'])
saveas(prec_vs_sharp,[FigPath 'prec_vs_sharp.pdf'])


