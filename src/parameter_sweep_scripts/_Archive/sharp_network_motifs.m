% script to search for and characterize sharp network motifs
clear 
close all
addpath(genpath('../utilities/'))

% define save path
ReadPath = ['../../out/bivariate_parameter_sweeps/'];
FigPath = '../../fig/sharpness_plots/';
mkdir(FigPath);

rate_bounds = [-4 -4 -4 -4 ;...
                4  4  4  4 ];%repmat([-4 ; 1],1,8); % constrain transition rate magnitude
              
rate_bounds = repmat(rate_bounds,1,2);              
[~,metric_names] = calculateMetrics_v4([]);

% get index of useful metrics
flux_index = find(strcmp(metric_names,'Flux'));
rate_index = find(strcmp(metric_names,'Production Rate'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
cycle_index = find(strcmp(metric_names,'CycleTime'));
state2_index = find(strcmp(metric_names,'2StateFrac'));
dKon_index = find(strcmp(metric_names,'dKondC'));
dKoff_index = find(strcmp(metric_names,'dKoffdC'));

% set basic plot parameters
markerSize = 50;
markerAlpha = 0.2;
pboc = [228,221,209]/255;
n_plot = 6e3;
% load simulation results
load_name_flux = ['param_sweep_results_' metric_names{flux_index} '_' ...
    metric_names{sharpness_index}];
  
load([ReadPath load_name_flux '_half_eq0.mat'])

%% perform clustering using normalize network rates
close all

network_rates = vertcat(sim_results_half.rate_array);
network_metrics = vertcat(sim_results_half.metric_array);
network_sharpness = network_metrics(:,sharpness_index);
network_cycle_times = network_metrics(:,cycle_index);

nanFilter = ~isnan(network_sharpness)err_filter = round(sim_struct_neq(2).metric_array(:,rate_index),2)~=0.50;;
network_rates = network_rates(nanFilter,:);
network_cycle_times = network_cycle_times(nanFilter,:);
network_sharpness = network_sharpness(nanFilter,:);

% set sharpness cutoff
sharpness_cutoff = prctile(network_sharpness,99);
sharpness_filter = network_sharpness>=sharpness_cutoff;

% get subset of sharp rates, normalized by cycle time 
network_rates_norm = network_rates./network_cycle_times;

% try simple PCA approach
[standard_data, mu, sigma] = zscore(network_rates_norm);     % standardize data so that the mean is 0 and the variance is 1 for each variable
[coeff, score, ~]  = pca(standard_data);     % perform PCA

plot_indices = randsample(1:size(score,1),n_plot,false);
% make figure
pca_fig = figure;

cmap = brewermap(9,'Set2');
hold on
s1 = scatter(score(plot_indices, 1), score(plot_indices, 2), markerSize ,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',...
                    markerAlpha, 'MarkerFaceColor',cmap(3,:));     
s1a = scatter(NaN,NaN,'MarkerEdgeAlpha',0, 'MarkerFaceColor', cmap(3,:));      
s2 = scatter(score(sharpness_filter, 1), score(sharpness_filter, 2), markerSize ,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',...
                    1, 'MarkerFaceColor',cmap(2,:)); 
               
box on    

legend([s1a, s2],'all networks','sharp networks','Location','southwest')

xlabel('1st principle component')
ylabel('2nd principle component')
set(gca,'FontSize',14)

ax = gca;

ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
set(gcf,'Color','w');
set(gca,'Color',pboc) 

pca_fig.InvertHardcopy = 'off';

saveas(pca_fig,[FigPath 'sharpness_pca.png'])
saveas(pca_fig,[FigPath 'sharpness_pca.pdf'])

%% Let's see if there is a clear trend in some other paramers
network_2state = network_metrics(nanFilter,state2_index);
network_dKon = network_metrics(nanFilter,dKon_index);
network_dKoff = network_metrics(nanFilter,dKoff_index);

% [KonEff,KoffEff,dKondC,dKoffdC] = calculateEffectiveRates24(network_rates,1);

weird_motifs = find(network_sharpness>0.495 & network_dKon > 1.9);
weird_networks = network_rates(weird_motifs,:);
weirdSS = fourStateOccupancy(weird_networks,1);
weirdSS(:,3).*weird_networks(:,1)
weirdSS(:,4).*weird_networks(:,5)


weirdSS(:,1).*weird_networks(:,3)
weirdSS(:,2).*weird_networks(:,5)

normal_motifs = find(network_sharpness>0.495 & network_dKon < 1.05);
normal_networks = network_rates(normal_motifs,:);
