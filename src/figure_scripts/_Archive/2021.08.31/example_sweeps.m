% script to generate data for system sharpness behaviors
clear 
close all
addpath(genpath('../utilities/'))
% define save path
OutPath = ['../../out/bivariate_parameter_sweeps/'];
mkdir(OutPath);

rate_bounds = [-4 -4 -4 -4 ;...
                4  4  4  4 ];%repmat([-4 ; 1],1,8); % constrain transition rate magnitude
              
rate_bounds = repmat(rate_bounds,1,2);              
[~,metric_names] = calculateMetrics_v4([]);

flux_index = find(strcmp(metric_names,'Flux')); % degree to which network breaks detailed balance
rate_index = find(strcmp(metric_names,'Production Rate')); % fraction of time spent in the active states (3 and 4)
sharpness_index = find(strcmp(metric_names,'Sharpness')); % derivative of production rate wrpt concncentration

% call edge sampler for production rate vs. sharpness
sample_indices = [sharpness_index,rate_index];


% set sim options
% n_seeds: number of mutants to generate from each boundary sample
% n_iters: number of resampling runs per simulation
% n_sim: number of unique edge sampling runs to do

sample_options = {'n_seeds',5,'n_iters',50,rate_bounds};

%% %%%%%%%%%%%%%%%% sharpness vs production rate %%%%%%%%%%%%%%%
% equilibrium_flag: if 1, reauire all networks to obey detailed balance

tic
sim_struct_neq = param_sweep_v5(sample_indices,sample_options{:},'equilibrium_flag',false);

sim_struct_eq = param_sweep_v5(sample_indices,sample_options{:},'equilibrium_flag',true);
toc

%% %%%%%%%%%%%%%%%% plot results %%%%%%%%%%%%%%%
markerAlpha = 0.3;
markerSize = 25;

% extract parameter vectors
metric_array_eq = sim_struct_eq.metric_array;        
mvec1_eq = reshape(metric_array_eq(:,sample_indices(1),:),1,[]);
mvec2_eq = reshape(metric_array_eq(:,sample_indices(2),:),1,[]);

metric_array_neq = sim_struct_neq.metric_array;        
mvec1_neq = reshape(metric_array_neq(:,sample_indices(1),:),1,[]);
mvec2_neq = reshape(metric_array_neq(:,sample_indices(2),:),1,[]);


example_figure = figure;

hold on    
cmap = brewermap(9,'Set2');

s2 = scatter(mvec1_neq,mvec2_neq,...
  markerSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(2,:));                

s1 = scatter(mvec1_eq,mvec2_eq,...
  markerSize,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(3,:));                

box on    

xlabel('sharpness');
ylabel('production rate');
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
ax.XAxis(1).FontSize = 11;


aff_rate.InvertHardcopy = 'off';
set(gcf,'color','w');

legend([s1 s2],'equilibrium','non-equilibrium','Color','w','Location','northwest')