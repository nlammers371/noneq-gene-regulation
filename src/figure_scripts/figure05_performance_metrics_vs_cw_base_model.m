% script to look at tradeoff between H0 and f

clear 
close all
addpath(genpath('../utilities/'))

% set basic parameters
nStates = 6;
% paramBounds = repmat([-10 ; 6],1,3*nStates-4); % constrain transition rate magnitude
[~,~,metric_names] = calculateMetricsSym_v2([]);
paramBounds = repmat([-5 ; 5],1,11); % constrain transition rate magnitude

% specify function path
% functionPath = ['../utilities/metricFunctions/n' num2str(nStates) '_OR/'];
functionPath = '../utilities/metricFunctions/symbolic/n006_s01_ns01_g01';

% make sure we're linked to the appropriate function subfolder% make sure we're linked to the appropriate function subfolder
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));

% define save path

DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\manuscript\';
% DropboxFolder = 'S:\Nick\Dropbox\Nonequilibrium\Nick\manuscript\';

FigPath = [DropboxFolder 'performance_metrics_vs_cw' filesep];
mkdir(FigPath);         

% get index of useful metrics
spec_index = find(strcmp(metric_names,'Specificity'));
% sharp_right_norm_index = find(strcmp(metric_names,'SharpnessRightNorm'));
ir_index = find(strcmp(metric_names,'DecisionRateNorm'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
precision_index = find(strcmp(metric_names,'Precision'));
precision_right_index = find(strcmp(metric_names,'PrecisionRight'));
sharp_right_alt_index = find(strcmp(metric_names,'SharpnessRight'));
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

cw = 1; % note that value does not matter
tic
[sim_info_neq, sim_struct_neq] = param_sweep_multi_v3([sharp_right_alt_index spec_index],functionPath,sweep_options{:},...
                                          'cw',cw,...
                                          'equilibrium_flag',false,'specFactor',alpha_factor,'paramBounds',paramBounds); 

[sim_info_eq, sim_struct_eq] = param_sweep_multi_v3([sharp_right_alt_index spec_index],functionPath,sweep_options{:},...
                                          'cw',cw,...
                                          'equilibrium_flag',true,'specFactor',alpha_factor,'paramBounds',paramBounds);                                   
toc     

% now for high CW
% cw_high = 1e3;
% tic
% [sim_info_neq_high, sim_struct_neq_high] = param_sweep_multi_v3([sharp_right_alt_index spec_index],functionPath,sweep_options{:},...
%                                           'half_max_flag',false,'cw',cw_high,...
%                                           'equilibrium_flag',false,'specFactor',alpha_factor);
% 
% [sim_info_eq_high, sim_struct_eq_high] = param_sweep_multi_v3([sharp_right_alt_index spec_index],functionPath,sweep_options{:},...
%                                           'half_max_flag',false,'cw',cw_high,...
%                                           'equilibrium_flag',true,'specFactor',alpha_factor);                                        
% toc     

%% Define function for analytic bound

% cr = 1;
f_axis = logspace(log10(alpha_factor),log10(alpha_factor^2));

% s0_bound2 = (1+(alpha_factor.^2+(-1).*f_axis).*((-1).*alpha_factor.*cw+(((-1)+alpha_factor).*cr+cw).*f_axis).^(-1));

close all
rd = [190 30 45]/255;

cr = 1;
rng(123)
cw_vec = [1 1e3];
high_flag = 0;
    
suffix = '';
if high_flag
    suffix = '_high';
end

cw = cw_vec(high_flag+1);
%     s0_bound_pd = (1 + (alpha_factor^2-f_axis)./(f_axis.*alpha_factor + cw*(f_axis-alpha_factor)));
s0_bound_pd = ((alpha_factor^2+alpha_factor.*f_axis-2*f_axis)./(f_axis.*alpha_factor - f_axis));
% generate vectors to plot (neq)
if ~high_flag
    results_array_neq = sim_struct_neq.metric_array;
else
    results_array_neq = sim_struct_neq_high.metric_array;
end
f0_scatter_vec_neq = 10.^results_array_neq(:,spec_index);
s0_scatter_vec_neq = results_array_neq(:,sharp_right_alt_index);
plot_options_neq = find(s0_scatter_vec_neq>=0 & f0_scatter_vec_neq >= 1);
plot_indices_neq = randsample(plot_options_neq,min([n_plot,length(plot_options_neq)]),false);

% generate vectors to plot (eq)
if ~high_flag
    results_array_eq = sim_struct_eq.metric_array;
else
    results_array_eq = sim_struct_eq_high.metric_array;
end
f0_scatter_vec_eq = 10.^results_array_eq(:,spec_index);
s0_scatter_vec_eq = results_array_eq(:,sharp_right_alt_index);
plot_options_eq = find(s0_scatter_vec_eq>=0 & f0_scatter_vec_eq >= 1);
plot_indices_eq = randsample(plot_options_eq,min([n_plot,length(plot_options_eq)]),false);

% make figure
s0_f0_fig = figure;
hold on
cmap = brewermap([],'Set2');

sneq = scatter(f0_scatter_vec_neq(plot_indices_neq)*alpha_factor,s0_scatter_vec_neq(plot_indices_neq),...
      markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha*0.5, 'MarkerFaceColor',cmap(2,:));

seq = scatter(f0_scatter_vec_eq(plot_indices_eq)*alpha_factor,s0_scatter_vec_eq(plot_indices_eq),...
      markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha*0.5, 'MarkerFaceColor',cmap(3,:));    

% plot predicted bound
plot(f_axis,s0_bound_pd,'-.','Color',rd,'LineWidth',3)

ylim([0 2])
xlim([alpha_factor alpha_factor^2])
set(gca,'xscale','log')    
xlabel('specificity (f)');
ylabel('intrinsic sharpness (S_0)')       
set(gca,'FontSize',14)
if ~high_flag
    set(gca,'Color',[228,221,209]/255) 
    grid on
end
set(gca,'xtick',[100 1e3 1e4]);%,'xticklabels',{'\alpha^0','\alpha^{0.5}','\alpha'})
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

s0_f0_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(s0_f0_fig,[FigPath 's0_vs_f0' suffix '.png'])
saveas(s0_f0_fig,[FigPath 's0_vs_f0' suffix '.pdf'])

