% Plot results for sharpness vs precision for higher order models
clear 
close all
addpath(genpath('../utilities/'))

% %%%%%%%%%%%%%%%%  set relative read and write paths %%%%%%%%%%%%%%%%%%%%
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\';
DataPath = [DropboxFolder  'SweepOutput\sweeps03_info_vs_cw' filesep ];
FigPath = [DropboxFolder '\manuscript\writeup' filesep];
mkdir(FigPath);

% %%%%%%%%%%%%%%%%  Set plot parameters and constants %%%%%%%%%%%%%%%%%%%%
n_plot = 3e3; % number of points to plot
n_samp = 5e4; % number of points to plot
markerAlpha = 0.5; % marker transparency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%% Compare "equilibrium" strategy of adding binding sites with 

% get metric names for numeric sweeps
[~,~,metric_names_num] = calculateMetricsNumeric_v3([]);

ir_index_num = find(strcmp(metric_names_num,'IR'));
cw_index_num = find(strcmp(metric_names_num,'CW'));
spec_index_num = find(strcmp(metric_names_num,'Specificity'));
inv_dt_index_num = find(strcmp(metric_names_num,'InverseDecisionTime'));

cw_index = linspace(-1,6,101);

%%%%%%%%%%%%%%%%%%%%%%%%5
%%% Load multi bs results
% get list of sweep results files with only 1 genera TF reaction
multi_bs_sweep_files_eq = dir([DataPath 'sweep_results*g01*eq*']);
multi_bs_info_files_eq = dir([DataPath 'sweep_info*g01*eq*']);

multi_bs_sweep_files_neq = dir([DataPath 'sweep_results*g01_cw1.mat']);
multi_bs_info_files_neq = dir([DataPath 'sweep_info*g01_cw1.mat']);

% load
master_struct_multi_bs = struct;
for f = 1:length(multi_bs_sweep_files_eq)
  
    % load eq files
    load([DataPath multi_bs_sweep_files_eq(f).name])
    load([DataPath multi_bs_info_files_eq(f).name])
    
    master_struct_multi_bs(f).sweep_results_eq = sim_results;
    master_struct_multi_bs(f).sweep_info_eq = sim_info;
        
end    

rng(321);


for i = 1:length(master_struct_multi_bs)

    metric_array = master_struct_multi_bs(i).sweep_results_eq.metric_array;       
        
    % extract vectors    
    cw_vec = metric_array(:,cw_index_num);
    ir_vec = metric_array(:,ir_index_num)*bit_factor;            
    dt_vec = 1./metric_array(:,inv_dt_index_num);            
    spec_vec = metric_array(:,spec_index_num);
 
    % identify boundary points
    ir_b_vec = NaN(1,length(cw_index)-1);
    dt_b_vec = NaN(1,length(cw_index)-1);   
    spec_b_vec = NaN(1,length(cw_index)-1);
    
    for c = 1:length(cw_index)-1
        cw_filter = find(cw_vec<cw_index(c+1)&cw_vec>=cw_index(c));
        [ir_b_vec(c), mi] = nanmax(ir_vec(cw_filter));
        dt_b_vec(c) = dt_vec(cw_filter(mi));    
        spec_b_vec(c) = spec_vec(cw_filter(mi));
    
    end

    % store    
    master_struct_multi_bs(i).cw_boundary_eq = 10.^cw_index(2:end);
    master_struct_multi_bs(i).ir_boundary_eq = imgaussfilt(ir_b_vec,2);   
    master_struct_multi_bs(i).dt_boundary_eq = imgaussfilt(dt_b_vec,2);   
    master_struct_multi_bs(i).spec_boundary_eq = imgaussfilt(spec_b_vec,2);       
end
  
%%%%%%%%%%%%%%%%%%%%%%%%5
% Load multi g results
% get list of sweep results files with only 1 genera TF reaction
multi_g_sweep_files = dir([DataPath 'sweep_results_s01_ns00_g0*_cw1.mat']);
multi_g_info_files = dir([DataPath 'sweep_info_s01_ns00_g0*_cw1.mat']);

% load
master_struct_multi_g = struct;
for f = 1:length(multi_g_sweep_files)
  
    load([DataPath multi_g_sweep_files(f).name])
    load([DataPath multi_g_info_files(f).name])
    
    master_struct_multi_g(f).sweep_results = sim_results;
    master_struct_multi_g(f).sweep_info = sim_info;
end    

rng(321);

for i = 1:length(master_struct_multi_g)
  
    % extract vectors
    metric_array = master_struct_multi_g(i).sweep_results.metric_array;   
    cw_vec = metric_array(:,cw_index_num);
    ir_vec = metric_array(:,ir_index_num)*bit_factor;    
    dt_vec = 1./metric_array(:,inv_dt_index_num); 
    spec_vec = metric_array(:,spec_index_num);
        
    % identify boundary points
    ir_b_vec = NaN(1,length(cw_index)-1);
    dt_b_vec = NaN(1,length(cw_index)-1);
    spec_b_vec = NaN(1,length(cw_index)-1);

    for c = 1:length(cw_index)-1
        cw_filter = find(cw_vec<cw_index(c+1)&cw_vec>=cw_index(c));
        [ir_b_vec(c), mi] = nanmax(ir_vec(cw_filter));
        dt_b_vec(c) = dt_vec(cw_filter(mi));      
        spec_b_vec(c) = spec_vec(cw_filter(mi));        
    end
    
    % store    
    master_struct_multi_g(i).cw_boundary = 10.^cw_index(2:end);
    master_struct_multi_g(i).ir_boundary = imgaussfilt(ir_b_vec,2);  
    master_struct_multi_g(i).dt_boundary = imgaussfilt(dt_b_vec,2);      
    master_struct_multi_g(i).spec_boundary = imgaussfilt(spec_b_vec,2);   
end

% load dataset with decision time ranges 
load('decision_limit_info.mat','decision_limit_info')

% generate decision time range vectors
dt_cw_vec = geomean(vertcat(decision_limit_info.cw_ub,decision_limit_info.cw_lb));

% combine mous and human
dt_cw_vec(end-1) = mean(dt_cw_vec(end-1:end));
dt_cw_vec = dt_cw_vec(1:end-1);
dt_cw_vec = dt_cw_vec([1:2 4]);

%% Plot sharpness, precision, and specificity gains vs. cw 
alpha_index = master_struct_multi_bs(1).sweep_info_eq.a_index;
alpha_factor = master_struct_multi_bs(1).sweep_results_eq.rate_array(1,alpha_index);

close all
markerSize = 50;

% set plot parameters
rng(231);
cmap2 = brewermap(8,'set2');
% Define colormaps for use throughout
cmap_pu = brewermap(8,'Purples');
cmap_rd = brewermap(8,'Reds');
cmap_bu = brewermap(8,'Blues');
cmap_gre = brewermap(8,'Greens');
cmap_gra = brewermap(8,'Greys');


color_ind = 5;
cmap_cmb = [cmap_pu(color_ind,:); cmap_gre(color_ind,:); cmap_rd(color_ind,:); ...
          cmap_bu(color_ind,:) ; cmap_gra(color_ind,:)];

      
% make figure
close all

f0_fig = figure;
hold on

slb = [];

% plot 5 BS neq system first
% cw_vec = master_struct_multi_bs(5).cw_boundary_eq;
% f0_vec_5bs = 10.^master_struct_multi_bs(5).spec_boundary_eq;
% 
% slb(1) = scatter(cw_vec,f0_vec_5bs*alpha_factor, markerSize,'o','MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',1,...
%           'MarkerFaceColor',cmap_cmb(5,:));

plot(logspace(0,5),repelem(alpha_factor,50),'-.','Color','k','LineWidth',1)

% 1 BS and multiple conformations     
for i = 1:length(master_struct_multi_g)-2
    cw_vec = master_struct_multi_g(i).cw_boundary;
    f0_vec = 10.^master_struct_multi_g(i).spec_boundary;
    
    slb(i) = scatter(cw_vec,f0_vec*alpha_factor, markerSize,'o','MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',1,...
                      'MarkerFaceColor',cmap_pu(2*i,:));
end

% plot reference lines
yl = get(gca,'ylim');
yl(end) = 1e8;
yl(1) = 1;
plot(repmat(dt_cw_vec,2,1),repmat(yl',1,3),'-','Color','k','LineWidth',1)

% slb = [slb(2:end) slb(1)];                
legend(slb,'N_a=1','N_a=2','N_a=3','Location','northwest')
xlabel('relative wrong factor concentration (w/c)');
ylabel('specificity (f)')
% grid on
% box on
set(gca,'FontSize',14)
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'xtick',[1  10 10^2 10^3 10^4 10^5])
set(gca,'ytick',[1  10^2  10^4 10^6 10^8],'yticklabels',{'\alpha^0','\alpha','\alpha^2','\alpha^3','\alpha^4'})
xlim([1 1e5])
ylim([1 1e8])

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
grid on
f0_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
% ylim([0.25 1.75])

saveas(f0_fig,[FigPath 'f0_vs_cw.png'])
saveas(f0_fig,[FigPath 'f0_vs_cw.pdf']) 
