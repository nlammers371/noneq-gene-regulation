% Plot results of sharpness parameter sweeps
% clear 
close all
addpath(genpath('../utilities/'))

% %%%%%%%%%%%%%%%%  set relative read and write paths %%%%%%%%%%%%%%%%%%%%
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\';
DataPath = [DropboxFolder  'SweepOutput\sweeps01_info_vs_energy' filesep ];
FigPath = [DropboxFolder '\manuscript\appendices' filesep 'sweep_algorithm' filesep];
mkdir(FigPath);

% %%%%%%%%%%%%%%%%  Set plot parameters and constants %%%%%%%%%%%%%%%%%%%%
n_plot = 3e3; % number of points to plot
markerAlpha = 0.5; % marker transparency
markerSize = 75; % marker size

% set sweep options
sweep_options = {'n_sim',1,'n_seeds',5,'n_iters_max',50};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%% Plot Phi versus vs IR for simple 4 state gene circuit %%%%%%%%%%%%
% set path 
nStates = 4;
functionPath = ['../../utilities/metricFunctions/n' num2str(nStates) '_OR/'];

% get metric names 
[~,~,metric_names_sym] = calculateMetricsSym_v2([]);
sharp_index = find(strcmp(metric_names_sym,'Sharpness'));
prec_index = find(strcmp(metric_names_sym,'Precision'));


% run symbolic sweep
tic
[sweep_info_sym4, sweep_results_sym4] = ...
                param_sweep_multi_v3([sharp_index prec_index],...
                                    functionPath,sweep_options{:},...
                                    'equilibrium_flag',false,...
                                    'TauCycleTime',1,'downsample_output',0); 
toc
%% calculate maximum info rates
metric_array_4 = sweep_results_sym4.metric_array;   
sharpness_vec = metric_array_4(:,sharp_index);
precision_vec = exp(metric_array_4(:,prec_index));

% make vector to indicate iter id
iter_id_vec = sweep_results_sym4.iter_id_vec;

% define subdir
sweep_fig_dir = [FigPath 'sweep_frames' filesep];
mkdir(sweep_fig_dir);

% make figures
% plot_indices = randsample(1:length(iter_id_vec),3e3,false);
% phi_vec_plot = phi_vec;%(plot_indices);
% info_vec_plot = info_vec;%(plot_indices);
iter_id_vec_plot = iter_id_vec;%(plot_indices);
close all
for n = sweep_info_sym4.n_iter_vec+1
  
    H_P_fig = figure;%('Visible','off');
    
    cmap = flipud(brewermap(sweep_info_sym4.n_iter_vec+1,'Spectral'));    
    hold on        
    scatter(sharpness_vec(iter_id_vec_plot<=n),precision_vec(iter_id_vec_plot<=n),...
          markerSize,cmap(iter_id_vec_plot(iter_id_vec_plot<=n),:),'filled','MarkerEdgeAlpha',.5,'MarkerEdgeColor','k','MarkerFaceAlpha',0.5) 
      

    xlabel('normalized sharpness (S)');
    ylabel('normalized precision (P)')
    grid on
    set(gca,'FontSize',14)
    set(gca,'Color',[228,221,209]/255) 

    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';

    xlim([0 2.1])
    ylim([0 1.1])
    H_P_fig.InvertHardcopy = 'off';
    set(gcf,'color','w');      
    saveas(H_P_fig,[sweep_fig_dir 'S_vs_P' sprintf('%03d',n) '.png'])    
    if n == sweep_info_sym4.n_iter_vec+1
        colormap(cmap);
        h = colorbar;
        ylabel(h, 'sweep step')
        saveas(H_P_fig,[sweep_fig_dir 'S_vs_P' sprintf('%03d',n) '_cb.png']) 
        saveas(H_P_fig,[sweep_fig_dir 'S_vs_P' sprintf('%03d',n) '_cb.pdf']) 
    end
end

%% Make area plot


area_fig = figure;

plot(sweep_results_sym4.areaVec,'Color','k','LineWidth',3)

xlabel('sweep iteration');
ylabel('sweep area metric')

grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';


area_fig.InvertHardcopy = 'off';
set(gcf,'color','w');      
saveas(area_fig,[FigPath 'area_plot_4state.png'])   
saveas(area_fig,[FigPath 'area_plot_4state.pdf'])   
