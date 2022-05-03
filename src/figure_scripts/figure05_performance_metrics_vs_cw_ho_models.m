% Plot results for sharpness vs precision for higher order models
clear 
close all
addpath(genpath('../utilities/'))

% %%%%%%%%%%%%%%%%  set relative read and write paths %%%%%%%%%%%%%%%%%%%%
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\';
DataPath = [DropboxFolder  'SweepOutput\sweeps03_info_vs_cw' filesep ];
FigPath = [DropboxFolder '\manuscript\performance_metrics_vs_cw' filesep];
mkdir(FigPath);

% %%%%%%%%%%%%%%%%  Set plot parameters and constants %%%%%%%%%%%%%%%%%%%%
n_plot = 3e3; % number of points to plot
n_samp = 5e4; % number of points to plot
markerAlpha = 0.5; % marker transparency
bit_factor = log2(exp(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%% Compare "equilibrium" strategy of adding binding sites with 

% get metric names for numeric sweeps
[~,~,metric_names_num] = calculateMetricsNumeric_v3([]);

ir_index_num = find(strcmp(metric_names_num,'IR'));
cw_index_num = find(strcmp(metric_names_num,'CW'));
rate_index_num = find(strcmp(metric_names_num,'ProductionRate'));
phi_index_num = find(strcmp(metric_names_num,'Phi'));
tau_index_num = find(strcmp(metric_names_num,'TauCycle'));
inv_dt_index_num = find(strcmp(metric_names_num,'InverseDecisionTime'));
sharp_right_index_num = find(strcmp(metric_names_num,'SharpnessRight'));
sharp_index_num = find(strcmp(metric_names_num,'Sharpness'));
spec_index_num = find(strcmp(metric_names_num,'Specificity'));
prec_index_num = find(strcmp(metric_names_num,'Precision'));

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
        
    % load neq files
    load([DataPath multi_bs_sweep_files_neq(f).name])
    load([DataPath multi_bs_info_files_neq(f).name])
    
    master_struct_multi_bs(f).sweep_results = sim_results;
    master_struct_multi_bs(f).sweep_info = sim_info;
end    

rng(321);


for i = 1:length(master_struct_multi_bs)
    for eq = 0:1
        if eq
            metric_array = master_struct_multi_bs(i).sweep_results_eq.metric_array;  
        else
            metric_array = master_struct_multi_bs(i).sweep_results.metric_array;  
        end
        
        % extract vectors    
        cw_vec = metric_array(:,cw_index_num);
        ir_vec = metric_array(:,ir_index_num)*bit_factor;            
        dt_vec = 1./metric_array(:,inv_dt_index_num);            
        spec_vec = metric_array(:,spec_index_num);
        sharp_vec = metric_array(:,sharp_index_num);
        sharp_right_vec = metric_array(:,sharp_right_index_num);
        prec_vec = metric_array(:,prec_index_num);
    
        % identify boundary points
        ir_b_vec = NaN(1,length(cw_index)-1);
        dt_b_vec = NaN(1,length(cw_index)-1);
        sharp_b_vec = NaN(1,length(cw_index)-1);
        sharp_right_b_vec = NaN(1,length(cw_index)-1);
        spec_b_vec = NaN(1,length(cw_index)-1);
        prec_b_vec = NaN(1,length(cw_index)-1);
        sharp_b_vec_max = NaN(1,length(cw_index)-1);
        prec_b_vec_max = NaN(1,length(cw_index)-1);
        for c = 1:length(cw_index)-1
            cw_filter = find(cw_vec<cw_index(c+1)&cw_vec>=cw_index(c));
            [ir_b_vec(c), mi] = nanmax(ir_vec(cw_filter));
            dt_b_vec(c) = dt_vec(cw_filter(mi));
            sharp_b_vec(c) = sharp_vec(cw_filter(mi));        
            sharp_right_b_vec(c) = sharp_right_vec(cw_filter(mi));
            spec_b_vec(c) = spec_vec(cw_filter(mi));
            prec_b_vec(c) = prec_vec(cw_filter(mi));
            sharp_b_vec_max(c) = nanmax(sharp_vec(cw_filter));
            prec_b_vec_max(c) = nanmax(prec_vec(cw_filter));
        end
        
        if eq
            % store    
            master_struct_multi_bs(i).cw_boundary_eq = 10.^cw_index(2:end);
            master_struct_multi_bs(i).ir_boundary_eq = imgaussfilt(ir_b_vec,2);   
            master_struct_multi_bs(i).dt_boundary_eq = imgaussfilt(dt_b_vec,2);   
            master_struct_multi_bs(i).sharp_boundary_eq = imgaussfilt(sharp_b_vec,2);   
            master_struct_multi_bs(i).sharp_right_boundary_eq = imgaussfilt(sharp_right_b_vec,2);   
            master_struct_multi_bs(i).spec_boundary_eq = imgaussfilt(spec_b_vec,2);   
            master_struct_multi_bs(i).prec_boundary_eq = imgaussfilt(prec_b_vec,2);  
            master_struct_multi_bs(i).sharp_max_boundary_eq = imgaussfilt(sharp_b_vec_max,2);   
            master_struct_multi_bs(i).prec_max_boundary_eq = imgaussfilt(prec_b_vec_max,2);  
        else
            master_struct_multi_bs(i).cw_boundary = 10.^cw_index(2:end);
            master_struct_multi_bs(i).ir_boundary = imgaussfilt(ir_b_vec,2);   
            master_struct_multi_bs(i).dt_boundary = imgaussfilt(dt_b_vec,2);   
            master_struct_multi_bs(i).sharp_boundary = imgaussfilt(sharp_b_vec,2);   
            master_struct_multi_bs(i).sharp_right_boundary = imgaussfilt(sharp_right_b_vec,2);   
            master_struct_multi_bs(i).spec_boundary = imgaussfilt(spec_b_vec,2);   
            master_struct_multi_bs(i).prec_boundary = imgaussfilt(prec_b_vec,2);   
        end
    end
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
    sharp_vec = metric_array(:,sharp_index_num);
    sharp_right_vec = metric_array(:,sharp_right_index_num);
    prec_vec = metric_array(:,prec_index_num);
        
    % identify boundary points
    ir_b_vec = NaN(1,length(cw_index)-1);
    dt_b_vec = NaN(1,length(cw_index)-1);
    sharp_b_vec = NaN(1,length(cw_index)-1);
    sharp_right_b_vec = NaN(1,length(cw_index)-1);
    spec_b_vec = NaN(1,length(cw_index)-1);
    prec_b_vec = NaN(1,length(cw_index)-1);
    for c = 1:length(cw_index)-1
        cw_filter = find(cw_vec<cw_index(c+1)&cw_vec>=cw_index(c));
        [ir_b_vec(c), mi] = nanmax(ir_vec(cw_filter));
        dt_b_vec(c) = dt_vec(cw_filter(mi));
        sharp_b_vec(c) = sharp_vec(cw_filter(mi));        
        sharp_right_b_vec(c) = sharp_right_vec(cw_filter(mi));
        spec_b_vec(c) = spec_vec(cw_filter(mi));
        prec_b_vec(c) = prec_vec(cw_filter(mi));
    end
    
    % store    
    master_struct_multi_g(i).cw_boundary = 10.^cw_index(2:end);
    master_struct_multi_g(i).ir_boundary = imgaussfilt(ir_b_vec,2);  
    master_struct_multi_g(i).dt_boundary = imgaussfilt(dt_b_vec,2);  
    master_struct_multi_g(i).sharp_boundary = imgaussfilt(sharp_b_vec,2);   
    master_struct_multi_g(i).sharp_right_boundary = imgaussfilt(sharp_right_b_vec,2);   
    master_struct_multi_g(i).spec_boundary = imgaussfilt(spec_b_vec,2);   
    master_struct_multi_g(i).prec_boundary = imgaussfilt(prec_b_vec,2); 
end

%%%%%%%%%%%%%%%%%%%%%%%%5
% Load 2BS multi g results

% get list of sweep results files with only 1 genera TF reaction
multi_g2_sweep_files = dir([DataPath 'sweep_results_s02_ns00_g0*_cw1.mat']);
multi_g2_info_files = dir([DataPath 'sweep_info_s02_ns00_g0*_cw1.mat']);

% load
master_struct_multi_g2 = struct;
for f = 1:length(multi_g2_sweep_files)
  
    load([DataPath multi_g2_sweep_files(f).name])
    load([DataPath multi_g2_info_files(f).name])
    
    master_struct_multi_g2(f).sweep_results = sim_results;
    master_struct_multi_g2(f).sweep_info = sim_info;
end    

rng(321);

for i = 1:length(master_struct_multi_g2)
  
    % extract vectors
    metric_array = master_struct_multi_g2(i).sweep_results.metric_array;   
    cw_vec = metric_array(:,cw_index_num);
    ir_vec = metric_array(:,ir_index_num)*bit_factor;    
    dt_vec = 1./metric_array(:,inv_dt_index_num); 
    spec_vec = metric_array(:,spec_index_num);
    sharp_vec = metric_array(:,sharp_index_num);
    sharp_right_vec = metric_array(:,sharp_right_index_num);
    prec_vec = metric_array(:,prec_index_num);
    
    % identify boundary points
    ir_b_vec = NaN(1,length(cw_index)-1);
    dt_b_vec = NaN(1,length(cw_index)-1);
    sharp_b_vec = NaN(1,length(cw_index)-1);
    sharp_right_b_vec = NaN(1,length(cw_index)-1);
    spec_b_vec = NaN(1,length(cw_index)-1);
    prec_b_vec = NaN(1,length(cw_index)-1);
    for c = 1:length(cw_index)-1
        cw_filter = find(cw_vec<cw_index(c+1)&cw_vec>=cw_index(c));
        [ir_b_vec(c), mi] = nanmax(ir_vec(cw_filter));
        dt_b_vec(c) = dt_vec(cw_filter(mi));
        sharp_b_vec(c) = sharp_vec(cw_filter(mi));        
        sharp_right_b_vec(c) = sharp_right_vec(cw_filter(mi));
        spec_b_vec(c) = spec_vec(cw_filter(mi));
        prec_b_vec(c) = prec_vec(cw_filter(mi));
    end
    
    % store    
    master_struct_multi_g2(i).cw_boundary = 10.^cw_index(2:end);
    master_struct_multi_g2(i).ir_boundary = imgaussfilt(ir_b_vec,2);  
    master_struct_multi_g2(i).dt_boundary = imgaussfilt(dt_b_vec,2);
    master_struct_multi_g2(i).sharp_boundary = imgaussfilt(sharp_b_vec,2);   
    master_struct_multi_g2(i).sharp_right_boundary = imgaussfilt(sharp_right_b_vec,2);   
    master_struct_multi_g2(i).spec_boundary = imgaussfilt(spec_b_vec,2);   
    master_struct_multi_g2(i).prec_boundary = imgaussfilt(prec_b_vec,2); 
end

%% load dataset with decision time ranges 
load('decision_limit_info.mat','decision_limit_info')

% generate decision time range vectors
dt_cw_vec = geomean(vertcat(decision_limit_info.cw_ub,decision_limit_info.cw_lb));

% combine mous and human
% dt_cw_vec(end-1) = mean(dt_cw_vec(end-1:end));
% dt_cw_vec = dt_cw_vec(1:end-1);
dt_cw_vec = dt_cw_vec([1:2 5]);

% Plot sharpness, precision, and specificity gains vs. cw 
close all
markerSize = 75;
ds_factor = 3;
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
cmap_cmb = [cmap_pu(3,:); cmap_gre(color_ind,:); cmap_rd(color_ind,:); ...
          cmap_bu(color_ind,:) ; cmap_gra(color_ind,:)];


red1 = [242 107 98]/255;
red2 = [250 187 161]/255;
red3 = mean(vertcat(red1,red2));

green1 = [118 195 118]/255;
green2 = [200 228 191]/255;
green3 = mean(vertcat(green1,green2));

s_p_fig = figure;
hold on

xlim([1 1e5])
ylim([1e-1 1e2])
% ylim([1e-6 1])

% plot reference lines
yl = get(gca,'ylim');
% yl(end) = 5;
p3 = plot(repmat(dt_cw_vec,2,1),repmat(yl',1,3),'-.','Color','k','LineWidth',1.5);

% 1 BS and multiple conformations     
for i = 1%:length(master_struct_multi_g)-1
    cw_vec = master_struct_multi_g(i).cw_boundary;
    ds_vec = 1:ds_factor:length(cw_vec);
    s_vec = master_struct_multi_g(i).sharp_boundary./master_struct_multi_bs(1).sharp_max_boundary_eq;
    p_vec = master_struct_multi_g(i).prec_boundary./master_struct_multi_bs(1).prec_max_boundary_eq;
    
    p1 = plot(cw_vec,s_vec, 'Color',green2,'LineWidth',2);
    s1 = scatter(cw_vec(ds_vec),s_vec(ds_vec), markerSize,'MarkerEdgeAlpha',1,'MarkerEdgeColor','k','MarkerFaceAlpha',1,...
                      'MarkerFaceColor',green3);
                    
    p2 = plot(cw_vec,p_vec, 'Color',red2,'LineWidth',2);
    s2 = scatter(cw_vec(ds_vec),p_vec(ds_vec), markerSize,'d','MarkerEdgeAlpha',1,'MarkerEdgeColor','k','MarkerFaceAlpha',1,...
                      'MarkerFaceColor',red3);
end
                    
% grid on
% box on
set(gca,'FontSize',14)
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'xtick',[1  10^1 10^2  10^3 10^4 10^5])


% slb = [slb(2:end) slb(1)];                
xlabel('relative wrong factor concentration (w/c)');
ylabel('nonequilibrium gain')

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
grid on
s_p_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
% ylim([0.25 1.75])

saveas(s_p_fig,[FigPath 's_p_vs_cw_LC2.png'])
saveas(s_p_fig,[FigPath 's_p_vs_cw_LC2.pdf']) 


delete([p1 p2 p3]);
delete([s1 s2]);

slb = [];
plb = [];

for i = 1:length(master_struct_multi_g)-1
    cw_vec = master_struct_multi_g(i).cw_boundary;
    ds_vec = 1:ds_factor:length(cw_vec);
    s_vec = master_struct_multi_g(i).sharp_boundary./master_struct_multi_bs(1).sharp_max_boundary_eq;
    p_vec = master_struct_multi_g(i).prec_boundary./master_struct_multi_bs(1).prec_max_boundary_eq;
    
    slb(end+1) = plot(cw_vec,s_vec, 'Color',cmap_pu(2*i,:),'LineWidth',2);
    scatter(cw_vec(ds_vec),s_vec(ds_vec), markerSize,'MarkerEdgeAlpha',1,'MarkerEdgeColor','k','MarkerFaceAlpha',1,...
                      'MarkerFaceColor',brighten(cmap_pu(2*i,:),0.1));
                    
    plot(cw_vec,p_vec, 'Color',cmap_pu(2*i,:),'LineWidth',2)
    scatter(cw_vec(ds_vec),p_vec(ds_vec), markerSize,'d','MarkerEdgeAlpha',1,'MarkerEdgeColor','k','MarkerFaceAlpha',1,...
                      'MarkerFaceColor',brighten(cmap_pu(2*i,:),0.1));
end
ylim([1e-1 1e3])
yl = get(gca,'ylim');
% yl(end) = 5;
p1 = plot(repmat(dt_cw_vec,2,1),repmat(yl',1,3),'-.','Color','k','LineWidth',1.5);
legend(slb, 'N_{LC}=2','N_{LC}=3','N_{LC}=4','N_{LC}=5','Location','northwest','Color','w')
saveas(s_p_fig,[FigPath 's_p_vs_cw.png'])
saveas(s_p_fig,[FigPath 's_p_vs_cw.pdf'])

%% get some numbers
i = 1;
cw_vec = master_struct_multi_g(i).cw_boundary;
ds_vec = 1:ds_factor:length(cw_vec);
s_vec = master_struct_multi_g(i).sharp_boundary./master_struct_multi_bs(1).sharp_max_boundary_eq;
p_vec = master_struct_multi_g(i).prec_boundary./master_struct_multi_bs(1).prec_max_boundary_eq;

p_gain = p_vec(find(cw_vec<=dt_cw_vec(end),1,'last'))
s_gain = s_vec(find(cw_vec<=dt_cw_vec(end),1,'last'))

%% s0 vs. f figure

close all

s0_f_fig = figure;
hold on

slb = [];

% 1 BS and multiple conformations     
for i = 1:length(master_struct_multi_g)-1
    cw_vec = master_struct_multi_g(i).cw_boundary;
    ds_vec = 1:ds_factor:length(cw_vec);
    s_vec = master_struct_multi_g(i).sharp_right_boundary;
    f0_vec = 10.^master_struct_multi_g(i).spec_boundary;
    
    slb(i) = plot(cw_vec,s_vec, '-','Color',cmap_pu(2*i,:),'LineWidth',2);
    scatter(cw_vec(ds_vec),s_vec(ds_vec), markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',1,...
                      'MarkerFaceColor',brighten(cmap_pu(2*i,:),0.1));
    
    yyaxis right
    plot(cw_vec,f0_vec,'-','Color',cmap_pu(2*i,:),'LineWidth',2)
    scatter(cw_vec(ds_vec),f0_vec(ds_vec), markerSize,'^','MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',1,...
                      'MarkerFaceColor',brighten(cmap_pu(2*i,:),0.1));
    
    yyaxis left
end
                    

% plot reference lines
ylim([10^(-2/8) 10^2])
yl = get(gca,'ylim');
set(gca,'yscale','log')
% yl(end) = 10;
plot(repmat(dt_cw_vec,2,1),repmat(yl',1,3),'-.','Color','k','LineWidth',1.5)
ylabel('intrinsic sharpness (S_0/N_b)')

% switch to f axis
yyaxis right
% grid on
ylim([0.1 1e8])
set(gca,'yscale','log')
ylabel('specificity (f/\alpha)')

% slb = [slb(2:end) slb(1)];                
legend(slb, 'N_{LC}=2','N_{LC}=2','N_{LC}=3','N_{LC}=4','Location','northwest')
xlabel('relative wrong factor concentration (w/c)');

set(gca,'FontSize',14)
set(gca,'xscale','log')
set(gca,'xtick',[1  10^2  10^4])
xlim([1 1e5])
% ylim([1e-6 1])


ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
ax.XAxis(1).Color = 'k';

ax.XGrid ='on';
ax.YGrid(1) ='on'; 

s0_f_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255)

saveas(s0_f_fig,[FigPath 's0_f_vs_cw.png'])
saveas(s0_f_fig,[FigPath 's0_f_vs_cw.pdf']) 

% 
% %%
% % make figure
% close all
% 
% s0_fig = figure;
% hold on
% 
% slb = [];
% 
% % plot 5 BS neq system first
% cw_vec = master_struct_multi_bs(5).cw_boundary;
% s_vec_5bs_eq = nanmax(master_struct_multi_bs(5).sharp_right_boundary_eq);
% s_vec_5bs = master_struct_multi_bs(5).sharp_right_boundary/s_vec_5bs_eq;
% 
% slb(1) = scatter(cw_vec,s_vec_5bs, markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',1,...
%           'MarkerFaceColor',cmap_cmb(5,:));
% 
% % 1 BS and multiple conformations     
% for i = 1:length(master_struct_multi_g)-1
%     cw_vec = master_struct_multi_g(i).cw_boundary;
%     s_vec = master_struct_multi_g(i).sharp_right_boundary;
%     
%     slb(i+1) = scatter(cw_vec,s_vec, markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',1,...
%                       'MarkerFaceColor',cmap_pu(2*i,:));
% end
% s_vec_5bs_eq = nanmax(master_struct_multi_bs(2).sharp_right_boundary_eq);
% cw_vec = master_struct_multi_g2(3).cw_boundary;
% s_vec = master_struct_multi_g2(3).sharp_right_boundary/s_vec_5bs_eq;
% 
% slb(end+1) = scatter(cw_vec,s_vec, markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',1,...
%                   'MarkerFaceColor',cmap_gre(4,:));
%                     
% 
% % plot reference lines
% yl = get(gca,'ylim');
% yl(end) = 5;
% plot(repmat(dt_cw_vec,2,1),repmat(yl',1,3),'-.','Color','k','LineWidth',1)
% 
% slb = [slb(2:end) slb(1)];                
% legend(slb, 'N_g=2 (N_{LC}=1)','N_g=3 (N_{LC}=1)','N_g=4 (N_{LC}=1)','N_g=5 (N_{LC}=1)','N_g=4 (N_{LC}=2)','N_g=2 (N_{LC}=5)')
% xlabel('relative wrong factor concentration (c_w/c^*)');
% ylabel('intrinsic sharpness gain (H_0/H_0^{eq})')
% % grid on
% % box on
% set(gca,'FontSize',14)
% set(gca,'xscale','log')
% set(gca,'xtick',[1  10^2  10^4])
% xlim([1 1e5])
% % ylim([1e-6 1])
% 
% ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.XAxis(1).Color = 'k';
% grid on
% s0_fig.InvertHardcopy = 'off';
% set(gcf,'color','w');
% set(gca,'FontSize',14)
% set(gca,'Color',[228,221,209]/255) 
% % ylim([0.25 1.75])
% 
% saveas(s0_fig,[FigPath 's0_vs_cw.png'])
% saveas(s0_fig,[FigPath 's0_vs_cw.pdf']) 
% 
% %%
% f0_fig = figure;
% hold on
% 
% slb = [];
% 
% % plot 5 BS neq system first
% cw_vec = master_struct_multi_bs(5).cw_boundary;
% f0_vec_5bs = 10.^master_struct_multi_bs(5).spec_boundary;
% 
% slb(1) = scatter(cw_vec,f0_vec_5bs, markerSize,'d','MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',1,...
%           'MarkerFaceColor',cmap_cmb(5,:));
% 
% % 1 BS and multiple conformations     
% for i = 1:length(master_struct_multi_g)-1
%     cw_vec = master_struct_multi_g(i).cw_boundary;
%     f0_vec = 10.^master_struct_multi_g(i).spec_boundary;
%     
%     slb(i+1) = scatter(cw_vec,f0_vec, markerSize,'d','MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',1,...
%                       'MarkerFaceColor',cmap_pu(2*i,:));
% end
% 
% cw_vec = master_struct_multi_g2(3).cw_boundary;
% f0_vec = 10.^master_struct_multi_g2(3).spec_boundary;
% 
% slb(end+1) = scatter(cw_vec,f0_vec, markerSize,'d','MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',1,...
%                   'MarkerFaceColor',cmap_gre(4,:));
%                     
% % plot reference lines
% yl = get(gca,'ylim');
% yl(end) = 1e7;
% yl(1) = 0.1;
% plot(repmat(dt_cw_vec,2,1),repmat(yl',1,3),'-.','Color','k','LineWidth',1)
% 
% % slb = [slb(2:end) slb(1)];                
% % legend(slb, 'N_g=2 (N_{LC}=1)','N_g=3 (N_{LC}=1)','N_g=4 (N_{LC}=1)','N_g=5 (N_{LC}=1)','N_g=4 (N_{LC}=2)','N_g=2 (N_{LC}=5)','Location','northwest')
% xlabel('relative wrong factor concentration (c_w/c^*)');
% ylabel('specificity gain (f/f^{eq})')
% % grid on
% % box on
% set(gca,'FontSize',14)
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% set(gca,'xtick',[1  10^2  10^4])
% set(gca,'ytick',[1  10^2  10^4 10^6],'yticklabels',{'\alpha^0','\alpha','\alpha^2','\alpha^3'})
% xlim([1 1e5])
% ylim([.1 1e7])
% 
% ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.XAxis(1).Color = 'k';
% grid on
% f0_fig.InvertHardcopy = 'off';
% set(gcf,'color','w');
% set(gca,'FontSize',14)
% set(gca,'Color',[228,221,209]/255) 
% % ylim([0.25 1.75])
% 
% saveas(f0_fig,[FigPath 'f0_vs_cw.png'])
% saveas(f0_fig,[FigPath 'f0_vs_cw.pdf']) 
% 
% %%
% p_fig = figure;
% hold on
% 
% slb = [];
% 
% % plot 5 BS neq system first
% cw_vec = master_struct_multi_bs(5).cw_boundary;
% p_vec_5bs = master_struct_multi_bs(5).prec_boundary / sqrt(1/2);
% 
% slb(1) = scatter(cw_vec,p_vec_5bs, markerSize,'^','MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',1,...
%           'MarkerFaceColor',cmap_cmb(5,:));
% 
% % 1 BS and multiple conformations     
% for i = 1:length(master_struct_multi_g)-1
%     cw_vec = master_struct_multi_g(i).cw_boundary;
%     p_vec = master_struct_multi_g(i).prec_boundary / sqrt(1/2);
%     
%     slb(i+1) = scatter(cw_vec,p_vec, markerSize,'^','MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',1,...
%                       'MarkerFaceColor',cmap_pu(2*i,:));
% end
% 
% cw_vec = master_struct_multi_g2(3).cw_boundary;
% p0_vec = master_struct_multi_g2(3).prec_boundary/sqrt(1/2);
% 
% slb(end+1) = scatter(cw_vec,p0_vec, markerSize,'^','MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',1,...
%                   'MarkerFaceColor',cmap_gre(4,:));
%                     
% 
% yl = get(gca,'ylim');
% yl(end) = 1.5;
% yl(1) = 0;
% plot(repmat(dt_cw_vec,2,1),repmat(yl',1,3),'-.','Color','k','LineWidth',1)
% slb = [slb(2:end) slb(1)]; 
% 
% % legend(slb, 'N_g=2 (N_{LC}=1)','N_g=3 (N_{LC}=1)','N_g=4 (N_{LC}=1)','N_g=5 (N_{LC}=1)','N_g=4 (N_{LC}=2)','N_g=2 (N_{LC}=5)','Location','northwest')
% xlabel('relative wrong factor concentration (c_w/c^*)');
% ylabel('precision gain (P_0/P_0^{eq})')
% % grid on
% % box on
% set(gca,'FontSize',14)
% set(gca,'xscale','log')
% % set(gca,'yscale','log')
% set(gca,'xtick',[1  10^2  10^4])
% xlim([1 1e5])
% % ylim([.1 1e7])
% 
% ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.XAxis(1).Color = 'k';
% grid on
% p_fig.InvertHardcopy = 'off';
% set(gcf,'color','w');
% set(gca,'FontSize',14)
% set(gca,'Color',[228,221,209]/255) 
% % ylim([0.25 1.75])
% 
% saveas(p_fig,[FigPath 'p_vs_cw.png'])
% saveas(p_fig,[FigPath 'p_vs_cw.pdf']) 
% 
% % %% Check to see if sharpness bound holds for larger systems
% % cr = master_struct_multi_bs(3).sweep_info.crs;
% % alpha_factor = master_struct_multi_bs(3).sweep_info.defaultValues(master_struct_multi_bs(3).sweep_info.a_index);
% % f0_vec = logspace(log10(alpha_factor),log10(alpha_factor^2));
% % 
% % cw_vec = logspace(0,5);
% % % define expression specific to 1bs case
% % s_bound_fun = @(f0) ((-1)+alpha_factor).^(-1).*(alpha_factor.^2+((-2)+alpha_factor).*f0).*(cw_vec+cr.*f0).^(-1);
% % % define a more general ex[ression
% % s_bound_gen_fun = @(f0,fm,seq) f0./(cw_vec/cr + f0) .* (1 + (fm-f0)./(fm).*(2*seq-1));
% % 
% % close all
% % 
% % s_bound_s01 = s_bound_fun(alpha_factor);
% % s_bound_s01_2 = s_bound_gen_fun(alpha_factor,alpha_factor^2,1);
% % s_bound_f01 = s_bound_fun(alpha_factor^2);
% % s_bound_f01_2 = s_bound_gen_fun(alpha_factor^2,alpha_factor^2,1);
% % 
% % base_case = figure;
% % hold on
% % scatter(10.^master_struct_multi_g(1).sweep_results.metric_array(:,cw_index_num),master_struct_multi_g(1).sweep_results.metric_array(:,sharp_index_num))
% % plot(cw_vec,s_bound_s01,'LineWidth',2)
% % plot(cw_vec,s_bound_s01_2,'--','LineWidth',2)
% % plot(cw_vec,s_bound_f01,'LineWidth',2)
% % plot(cw_vec,s_bound_f01_2,'--','LineWidth',2)
% % set(gca,'xscale','log')
% % xlim([1 1e5])
% % 
% % %%
% % % let's see if something similar is at play in larger systems
% % cw_vec3 = master_struct_multi_bs(3).sweep_results.metric_array(:,cw_index_num);
% % f0_max3 = 100*10.^nanmax(master_struct_multi_bs(3).sweep_results.metric_array(cw_vec3>=3,spec_index_num));
% % seq_max3 = nanmax(master_struct_multi_bs(3).sweep_results_eq.metric_array(:,sharp_right_index_num));
% % 
% % s_bound_s03 = s_bound_gen_fun(alpha_factor,f0_max3,seq_max3);
% % s_bound_f03 = s_bound_gen_fun(f0_max3,f0_max3,seq_max3);
% % 
% % case3 = figure;
% % hold on
% % scatter(10.^master_struct_multi_bs(3).sweep_results.metric_array(:,cw_index_num),master_struct_multi_bs(3).sweep_results.metric_array(:,sharp_index_num))
% % plot(cw_vec,s_bound_s03,'LineWidth',2)
% % plot(cw_vec,s_bound_f03,'LineWidth',2)
% % set(gca,'xscale','log')
% % xlim([1 1e5])
% % 
% % %%
% % % let's see if something similar is at play in larger systems
% % cw_vec3G = master_struct_multi_g(3).sweep_results.metric_array(:,cw_index_num);
% % f0_max3G = 100*10.^nanmax(master_struct_multi_g(3).sweep_results.metric_array(cw_vec3>=3,spec_index_num));
% % seq_max3G = 1;%nanmax(master_struct_multi_g(3).sweep_results_eq.metric_array(:,sharp_right_index_num));
% % 
% % s_bound_s03G = 3/2*s_bound_gen_fun(alpha_factor,f0_max3G,seq_max3G);
% % s_bound_f03G = s_bound_gen_fun(f0_max3G,f0_max3G,seq_max3G);
% % 
% % case3 = figure;
% % hold on
% % scatter(10.^master_struct_multi_g(3).sweep_results.metric_array(:,cw_index_num),master_struct_multi_g(3).sweep_results.metric_array(:,sharp_index_num))
% % plot(cw_vec,s_bound_s03G,'LineWidth',2)
% % plot(cw_vec,s_bound_f03G,'LineWidth',2)
% % set(gca,'xscale','log')
% % xlim([1 1e5])