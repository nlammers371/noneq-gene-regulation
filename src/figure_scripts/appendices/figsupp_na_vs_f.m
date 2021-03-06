% Plot results for sharpness vs precision for higher order models
clear 
close all
addpath(genpath('../utilities/'))

% %%%%%%%%%%%%%%%%  set relative read and write paths %%%%%%%%%%%%%%%%%%%%
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\';
DataPath = [DropboxFolder  'SweepOutput\sweeps04B_f_vs_r_prec12' filesep ];
FigPath = [DropboxFolder '\manuscript\appendices' filesep];
mkdir(FigPath);

% %%%%%%%%%%%%%%%%  Set plot parameters and constants %%%%%%%%%%%%%%%%%%%%
n_samp = 5e4; % number of points to plot
markerAlpha = 0.5; % marker transparency
markerSize = 100; % marker size


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%% Look at impact of multiple molecular steps %%%%%%%%%%%%%%%%
[~,~,metric_names_num] = calculateMetricsNumeric_v3([]);
spec_index = find(strcmp(metric_names_num,'Specificity'));

% get list of sweep results files with only 1 genera TF reaction
multi_g_sweep_files = dir([DataPath 'sweep_results_s01_ns00_g0*neq*']);
multi_g_info_files = dir([DataPath 'sweep_info_s01_ns00_g0*neq*']);
  
% load
master_struct_multi_g = struct;
for f = 1:length(multi_g_sweep_files)
  
    load([DataPath multi_g_sweep_files(f).name])
    load([DataPath multi_g_info_files(f).name])
    
    master_struct_multi_g(f).sweep_results = sim_results;
    master_struct_multi_g(f).sweep_info = sim_info;
end    

alpha_factor = sim_info.specFactor;
rng(321);

xlg = [0.5 4.5];
close all

cmap_pu = brewermap(8,'Purples');
spec_a = figure;

hold on

% make plots
n_plot = 500;
x_vec = linspace(xlg(1), xlg(2));
plot(x_vec,alpha_factor.^(1+x_vec),':','Color','k','LineWidth',2)
for i = 1:length(master_struct_multi_g)
    f_vec = alpha_factor*10.^master_struct_multi_g(i).sweep_results.metric_array(:,spec_index);
    s_options = 1:length(f_vec);
    plot_ids = randsample(s_options,n_plot,true);
    g_temp = i + normrnd(0,0.025,n_plot,1);
    scatter(g_temp,f_vec(s_options(plot_ids)),25,'MarkerFaceColor',cmap_pu(i*2,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.1);
end

ylabel('non-equilibrium specificity (f)')
% ylim([2 5])
xlim(xlg)
xlabel('number of activation steps (N_{A})');
set(gca,'yscale','log')
grid on
box on
set(gca,'xtick',1:4)
set(gca,'ytick',[1 alpha_factor alpha_factor^2 alpha_factor^3 alpha_factor^4 alpha_factor^5],'yticklabel',{'\alpha^0','\alpha^1','\alpha^2','\alpha^3','\alpha^4','\alpha^5'} )
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
% ylim([0 0.4])
spec_a.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(spec_a,[FigPath 'spec_vs_a_sc.png'])
saveas(spec_a,[FigPath 'spec_vs_a_sc.pdf']) 



