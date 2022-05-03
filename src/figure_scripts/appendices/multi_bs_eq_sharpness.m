
clear 
close all
addpath(genpath('../utilities/'))

% %%%%%%%%%%%%%%%%  set relative read and write paths %%%%%%%%%%%%%%%%%%%%
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\';
DataPath = [DropboxFolder  'SweepOutput\sweeps02_sharpness_vs_precision_v2' filesep ];
DataPathSupp = [DropboxFolder  'SweepOutput\sweeps01_info_vs_energy_v3' filesep ];
FigPath = [DropboxFolder '\manuscript\Appendices' filesep];
mkdir(FigPath);

% %%%%%%%%%%%%%%%%  Set plot parameters and constants %%%%%%%%%%%%%%%%%%%%
n_plot = 3e3; % number of points to plot
n_samp = 5e4; % number of points to plot
markerAlpha = 0.5; % marker transparency
markerSize = 100; % marker size

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%% Examine impact of adding binding sites %%%%%%%%%%%%%%%%%
% get metric names for numeric sweeps
[~,~,metric_names_num] = calculateMetricsNumeric_v3([]);

ir_index_num = find(strcmp(metric_names_num,'IR'));
sharpness_index_num = find(strcmp(metric_names_num,'Sharpness'));
precision_index_num = find(strcmp(metric_names_num,'Precision'));
rate_index_num = find(strcmp(metric_names_num,'ProductionRate'));
phi_index_num = find(strcmp(metric_names_num,'Phi'));
tau_index_num = find(strcmp(metric_names_num,'TauCycle'));


% get list of sweep results files with only 1 genera TF reaction
multi_bs_sweep_files = dir([DataPath 'sweep_results*g01_cw0_eq*']);
multi_bs_info_files = dir([DataPath 'sweep_info*g01_cw0_eq*']);

% load
master_struct_multi_bs = struct;
for f = 1:length(multi_bs_sweep_files)
  
    load([DataPath multi_bs_sweep_files(f).name])
    load([DataPath multi_bs_info_files(f).name])
    
    master_struct_multi_bs(f).sweep_results = sim_results;
    master_struct_multi_bs(f).sweep_info = sim_info;
end    

% calculate sharpness vs precision boundaries 
% use precision and sharpness to eliminate outliers due to numerical 
% precision for now. Need to revisit this
p_max = [2.25 2.25 2.25 2.25 1.62.^2];
s_max = 100;

rng(321);
n_ir = 100;

for i = 1:length(master_struct_multi_bs)
    % extract vectors
    metric_array = master_struct_multi_bs(i).sweep_results.metric_array;   
    sharpness_vec = metric_array(:,sharpness_index_num).^2;
    precision_vec = metric_array(:,precision_index_num).^2;
    ir_vec = metric_array(:,ir_index_num);
    bound_filter = sharpness_vec<=s_max&precision_vec<=p_max(i);
    if i > 3
        ft2 = sqrt(sharpness_vec)<=0.15 & precision_vec >= 1;
    else
        ft2 = false(size(sharpness_vec));
    end
    options = find(bound_filter&~ft2);
    
    % identify boundary points
    b_points = boundary(sqrt(sharpness_vec(options)),sqrt(precision_vec(options)),0.9);
    
    % store
    master_struct_multi_bs(i).bound_points = options(b_points);
    master_struct_multi_bs(i).sharpness_boundary = sharpness_vec(options(b_points));
    master_struct_multi_bs(i).precision_boundary = precision_vec(options(b_points));
    
    % select N top performers to plot
    ir_99 = 0.99*nanmax(ir_vec(options));
    ir_options = find(ir_vec(options)>=ir_99);
    ir_indices = randsample(options(ir_options),n_ir,true);
    master_struct_multi_bs(i).sharp_ir = sharpness_vec(ir_indices);
    master_struct_multi_bs(i).prec_ir = precision_vec(ir_indices);
    
    % save max sharpness and precision
    master_struct_multi_bs(i).sharp_max = nanmax(sharpness_vec(options));
    master_struct_multi_bs(i).prec_max = nanmax(precision_vec(options));
end

%% Plot S vs. n_b

% Define colormaps for use throughout
cmap_pu = brewermap(8,'Purples');
cmap_rd = brewermap(8,'Reds');
cmap_bu = brewermap(8,'Blues');
cmap_gre = brewermap(8,'Greens');
cmap_gra = brewermap(8,'Greys');

close all
color_ind = 5;
cmap = [cmap_pu(3,:); cmap_gre(color_ind,:); cmap_rd(color_ind,:); ...
          cmap_bu(color_ind,:) ; cmap_gra(color_ind,:)];
       
xlg = [0.5 5.5];

close all
sharp_g = figure;

hold on

% make plots
n_plot = 500;
plot(linspace(xlg(1), xlg(2)),linspace(xlg(1), xlg(2)),':','Color','k','LineWidth',2)
for i = 1:length(master_struct_multi_bs)
    s_vec = master_struct_multi_bs(i).sweep_results.metric_array(:,sharpness_index_num);
    s_options = find(s_vec>=0);
    plot_ids = randsample(s_options,n_plot,true,ones(size(sqrt(s_vec(s_options)))));
    g_temp = i + normrnd(0,0.025,n_plot,1);
    scatter(g_temp,s_vec(plot_ids),25,'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.1);
end

ylabel('equilibrium sharpness (S)')
% ylim([2 5])
xlim(xlg)
xlabel('number of binding sites');

grid on
box on
set(gca,'xtick',1:4)
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
% ylim([0 0.4])
sharp_g.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(sharp_g,[FigPath 'sharp_vs_bs_eq.png'])
saveas(sharp_g,[FigPath 'sharp_vs_bs_eq.pdf']) 
