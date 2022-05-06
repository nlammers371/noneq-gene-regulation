% Plot results for IR vs energy for higher order models
clear 
close all
addpath(genpath('../utilities/'))

% %%%%%%%%%%%%%%%%  set relative read and write paths %%%%%%%%%%%%%%%%%%%%
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\';
DataPath_a = [DropboxFolder  'SweepOutput\sweeps01_info_vs_energy_comparison' filesep ];
DataPath_b = [DropboxFolder  'SweepOutput\sweeps01_info_vs_energy' filesep ];

FigPath = [DropboxFolder '\manuscript\appendices' filesep 'sweep_algorithm' filesep];
mkdir(FigPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% load numerical sweep results sets MULTI BS

%%%%%%%%%
% Multi BS
%%%%%%%%%

% get list of sweep results files with only 1 genera TF reaction
multi_bs_sweep_files_a = dir([DataPath_a 'sweep_results*g01*']);
multi_bs_info_files_a = dir([DataPath_a 'sweep_info*g01*']);

multi_bs_sweep_files_b = dir([DataPath_b 'sweep_results*g01*']);
multi_bs_info_files_b = dir([DataPath_b 'sweep_info*g01*']);

% load
master_struct_multi_bs = struct;
for f = 1:length(multi_bs_info_files_a)
  
    % load "a" files
    load([DataPath_a  multi_bs_info_files_a(f).name])
    load([DataPath_a  multi_bs_sweep_files_a(f).name])
    
    master_struct_multi_bs(f).sweep_results_a = sim_results;
    master_struct_multi_bs(f).sweep_info_a = sim_info;
    
    clear sim_results
    clear sim_info
    
    % load "b" files
    load([DataPath_b  multi_bs_info_files_b(f).name])
    load([DataPath_b  multi_bs_sweep_files_b(f).name])
    
    master_struct_multi_bs(f).sweep_results_b = sim_results;
    master_struct_multi_bs(f).sweep_info_b = sim_info;
    
    clear sim_results
    clear sim_info
    
end            



%%%%%%%%%
% Multi-LC
%%%%%%%%%

% get list of sweep results files with only 1 genera TF reaction
multi_lc_sweep_files_a = dir([DataPath_a 'sweep_results*s01*']);
multi_lc_info_files_a = dir([DataPath_a 'sweep_info*s01*']);

multi_lc_sweep_files_b = dir([DataPath_b 'sweep_results*s01*']);
multi_lc_info_files_b = dir([DataPath_b 'sweep_info*s01*']);

% load
master_struct_multi_lc = struct;
for f = 1:4
  
    load([DataPath_a  multi_lc_sweep_files_a(f).name])
    load([DataPath_a  multi_lc_info_files_a(f).name])
    
    master_struct_multi_lc(f).sweep_results_a = sim_results;
    master_struct_multi_lc(f).sweep_info_a = sim_info;
    
    clear sim_results
    clear sim_info
    
    load([DataPath_b  multi_lc_sweep_files_b(f).name])
    load([DataPath_b  multi_lc_info_files_b(f).name])
    
    master_struct_multi_lc(f).sweep_results_b = sim_results;
    master_struct_multi_lc(f).sweep_info_b = sim_info;
    
    clear sim_results
    clear sim_info
end    

%% Calculate upper information bounds as a function of Phi for all architectures
minDP = 100;
[~,~,metric_names_num] = calculateMetricsNumeric_v3([]);
ir_index = find(strcmp(metric_names_num,'IR'));
phi_index = find(strcmp(metric_names_num,'Phi'));
n_groups = 31;
phi_bins = logspace(-1,log10(5000),n_groups+1);

% multi BS
ir_max_array_bs = NaN(length(phi_bins),length(master_struct_multi_bs),2);

for i = 1:length(master_struct_multi_bs)
    % extract metric vectors
    metric_array_a = master_struct_multi_bs(i).sweep_results_a.metric_array;
    metric_array_b = master_struct_multi_bs(i).sweep_results_b.metric_array;
    
    ir_vec_a = metric_array_a(:,ir_index)*log2(exp(1));
    phi_vec_a = metric_array_a(:,phi_index);
    
    ir_vec_b = metric_array_b(:,ir_index)*log2(exp(1));
    phi_vec_b = metric_array_b(:,phi_index);
    
    for p = 1:n_groups
        phi_filter_a = phi_vec_a>=phi_bins(p) & phi_vec_a<phi_bins(p+1);
        phi_filter_b = phi_vec_b>=phi_bins(p) & phi_vec_b<phi_bins(p+1);
        
        if sum(phi_filter_a)>=minDP && sum(phi_filter_b)>=minDP
            ir_max_array_bs(p,i,1) = nanmax(ir_vec_a(phi_filter_a));
            ir_max_array_bs(p,i,2) = nanmax(ir_vec_b(phi_filter_b));
        end
    end
end    

% multi LC
ir_max_array_lc = NaN(length(phi_bins),length(master_struct_multi_lc),2);

for i = 1:length(master_struct_multi_lc)
    % extract metric vectors
    metric_array_a = master_struct_multi_lc(i).sweep_results_a.metric_array;
    metric_array_b = master_struct_multi_lc(i).sweep_results_b.metric_array;
    
    ir_vec_a = metric_array_a(:,ir_index)*log2(exp(1));
    phi_vec_a = metric_array_a(:,phi_index);
    
    ir_vec_b = metric_array_b(:,ir_index)*log2(exp(1));
    phi_vec_b = metric_array_b(:,phi_index);
    
    for p = 1:n_groups
        phi_filter_a = phi_vec_a>=phi_bins(p) & phi_vec_a<phi_bins(p+1);
        phi_filter_b = phi_vec_b>=phi_bins(p) & phi_vec_b<phi_bins(p+1);
        
        if sum(phi_filter_a)>=minDP && sum(phi_filter_b)>=minDP
            ir_max_array_lc(p,i,1) = nanmax(ir_vec_a(phi_filter_a));
            ir_max_array_lc(p,i,2) = nanmax(ir_vec_b(phi_filter_b));
        end
    end
end    
    
%%
ir_match_lc = figure;
hold on
cmap_pu = brewermap(8,'Purples');

s = [];
for i = length(master_struct_multi_lc):-1:1
    s(end+1) = scatter(ir_max_array_lc(:,i,1),ir_max_array_lc(:,i,2),75,'MarkerFaceColor',cmap_pu(2+i,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',0.75);
end

xlabel('maximum IR value (sweep round a)');
ylabel('maximum IR value (sweep round b)');

grid on
set(gca,'FontSize',14)
legend(fliplr(s),'N_{LC}=2','N_{LC}=3','N_{LC}=4','N_{LC}=5','Location','southeast')
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';


ir_match_lc.InvertHardcopy = 'off';
set(gcf,'color','w');      
saveas(ir_match_lc,[FigPath 'ir_phi_scatter_LC.png'])   
saveas(ir_match_lc,[FigPath 'ir_phi_scatter_LC.pdf']) 

%% 
% Define colormaps for use throughout
cmap_pu = brewermap(8,'Purples');
cmap_rd = brewermap(8,'Reds');
cmap_bu = brewermap(8,'Blues');
cmap_gre = brewermap(8,'Greens');
cmap_gra = brewermap(8,'Greys');

close all
color_ind = 5;
cmap_full = [cmap_pu(3,:); cmap_gre(color_ind,:); cmap_rd(color_ind,:); ...
          cmap_bu(color_ind,:) ; cmap_gra(color_ind,:)];

        
ir_match_bs = figure;
hold on
cmap_pu = brewermap(8,'Purples');

s = [];
for i = length(master_struct_multi_bs):-1:1
    s(end+1) = scatter(ir_max_array_bs(:,i,1),ir_max_array_bs(:,i,2),75,'MarkerFaceColor',cmap_full(i,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',0.75);
end

xlabel('maximum IR value (sweep round a)');
ylabel('maximum IR value (sweep round b)');

grid on
set(gca,'FontSize',14)
legend(fliplr(s),'N_{B}=1','N_{B}=2','N_{B}=3','N_{B}=4','N_{B}=5','Location','southeast')
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';


ir_match_bs.InvertHardcopy = 'off';
set(gcf,'color','w');      
saveas(ir_match_bs,[FigPath 'ir_phi_scatter_BS.png'])   
saveas(ir_match_bs,[FigPath 'ir_phi_scatter_BS.pdf']) 