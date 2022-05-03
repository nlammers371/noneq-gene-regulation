% Plot results of sharpness parameter sweeps
% clear 
close all
addpath(genpath('../utilities/'))

% %%%%%%%%%%%%%%%%  set relative read and write paths %%%%%%%%%%%%%%%%%%%%
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\';
DataPath = [DropboxFolder  'SweepOutput\sweeps01_info_vs_energy' filesep ];
FigPath = [DropboxFolder '\manuscript\appendices' filesep];
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
functionPath = ['../utilities/metricFunctions/n' num2str(nStates) '_OR/'];

% get metric names 
[~,~,metric_names_sym] = calculateMetricsSym_v2([]);
ir_index = find(strcmp(metric_names_sym,'DecisionRateNorm'));
decision_time_index = find(strcmp(metric_names_sym,'DecisionTimeNorm'));
sharpness_index = find(strcmp(metric_names_sym,'Sharpness'));
cycle_time_index = find(strcmp(metric_names_sym,'CycleTime'));
phi_index = find(strcmp(metric_names_sym,'Phi'));

% run symbolic sweep
tic
[sweep_info_sym4, sweep_results_sym4] = ...
                param_sweep_multi_v3([phi_index ir_index],...
                                    functionPath,sweep_options{:},...
                                    'half_max_flag',false,'equilibrium_flag',false,...
                                    'TauCycleTime',1,'downsample_output',0); 
toc
%% calculate maximum info rates
metric_array_phi = vertcat(sweep_results_sym4.metric_array);   
phi_vec = metric_array_phi(:,phi_index);
sharpness_vec = metric_array_phi(:,sharpness_index);
info_vec = metric_array_phi(:,ir_index) * log2(exp(1));% ./ c_factor;


% make vector to indicate iter id
iter_id_vec = sweep_results_sym4.iter_id_vec;

% define subdir
sweep_fig_dir = [FigPath 'sweep_frames' filesep];
mkdir(sweep_fig_dir);
% make figures
plot_indices = randsample(1:length(iter_id_vec),3e3,false);
phi_vec_plot = phi_vec;%(plot_indices);
info_vec_plot = info_vec;%(plot_indices);
iter_id_vec_plot = iter_id_vec;%(plot_indices);
close all
for n = sweep_info_sym4.n_iter_vec+1
    info_vs_flux = figure;
    if n <= sweep_info_sym4.n_iter_vec
        cmap = flipud(brewermap(sweep_info_sym4.n_iter_vec,'Spectral'));
        hold on        
        scatter(phi_vec_plot(iter_id_vec_plot<=n)/5,info_vec_plot(iter_id_vec_plot<=n)/5,...
              markerSize,cmap(iter_id_vec_plot(iter_id_vec_plot<=n),:),'filled','MarkerEdgeAlpha',.5,'MarkerEdgeColor','k','MarkerFaceAlpha',0.5) 
    else
      cmap = brewermap(8,'Purples');
      hold on
      
      scatter(phi_vec/5,info_vec/5,...
            markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',0.75,...
            'MarkerFaceColor',cmap(5,:)); 
          
      phi_axis = logspace(-3,log10(40));
      ir_max_vec = NaN(1,length(phi_axis)-1);
      for p = 1:length(phi_axis)-1
          ir_max_vec(p) = nanmax(info_vec(phi_vec<phi_axis(p+1)&phi_vec>=phi_axis(p)))/5;
      end
      plot(phi_axis(2:end)/5,imgaussfilt(ir_max_vec,1.5),'Color','k','LineWidth',3)
      
    end
    ylim([0 3e-3])
    xlim([5e-3 3e1]/5)

    xlabel('energy dissipation rate (k_BT per minute)');
    ylabel('information rate (bits per minute)')
    grid on
    set(gca,'FontSize',14)
    set(gca,'Color',[228,221,209]/255) 

    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';

    info_vs_flux.InvertHardcopy = 'off';
    set(gcf,'color','w');
    set(gca,'xscale','log')
    if n <= sweep_info_sym4.n_iter_vec
        saveas(info_vs_flux,[sweep_fig_dir 'info_vs_flux' sprintf('%03d',n) '.tif'])    
    else
        saveas(info_vs_flux,[FigPath 'info_vs_flux.tif'])    
    end
end
