% Plot results for sharpness vs precision for higher order models
clear 
close all
addpath(genpath('../utilities/'))

% %%%%%%%%%%%%%%%%  set relative read and write paths %%%%%%%%%%%%%%%%%%%%
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\';
DataPath = [DropboxFolder  'SweepOutput\sweeps04_s0_vs_f0' filesep ];
FigPath = [DropboxFolder '\manuscript\performance_metrics_vs_cw' filesep];
mkdir(FigPath);

% %%%%%%%%%%%%%%%%  Set plot parameters and constants %%%%%%%%%%%%%%%%%%%%
n_plot = 3e3; % number of points to plot
n_samp = 5e4; % number of points to plot
markerAlpha = 0.5; % marker transparency
markerSize = 75; % marker size
bit_factor = log2(exp(1));
alpha_factor = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%% Compare "equilibrium" strategy of adding binding sites with 

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

spec_index = linspace(-0.1,8,201);

%%%%%%%%%%%%%%%%%%%%%%%%5
%%% Load multi bs results
% get list of sweep results files with only 1 genera TF reactio
multi_bs_sweep_files_neq = dir([DataPath 'sweep_results*g01_cw1_neq.mat']);
multi_bs_info_files_neq = dir([DataPath 'sweep_info*g01_cw1_neq.mat']);

multi_bs_sweep_files_neq_high = dir([DataPath 'sweep_results*g01_cw1_neq_cval1000.mat']);
multi_bs_info_files_neq_high = dir([DataPath 'sweep_info*g01_cw1_neq_cval1000.mat']);

% load
master_struct_multi_bs = struct;
for f = 1:length(multi_bs_sweep_files_neq)
        
    % load neq files (cw=1)
    load([DataPath multi_bs_sweep_files_neq(f).name])
    load([DataPath multi_bs_info_files_neq(f).name])
    
    master_struct_multi_bs(f).sweep_results = sim_results;
    master_struct_multi_bs(f).sweep_info = sim_info;
    
    % load neq files (cw=1)
    load([DataPath multi_bs_sweep_files_neq_high(f).name])
    load([DataPath multi_bs_sweep_files_neq_high(f).name])
    
    master_struct_multi_bs(f).sweep_results_high = sim_results;
    master_struct_multi_bs(f).sweep_info_high = sim_info;
end    

rng(321);


for i = 1:length(master_struct_multi_bs)
    for high_flag = 0:1
        if ~high_flag
            metric_array = master_struct_multi_bs(i).sweep_results.metric_array;          
        else
            metric_array = master_struct_multi_bs(i).sweep_results_high.metric_array;          
        end
        % extract vectors        
        spec_vec = metric_array(:,spec_index_num);    
        prec_vec = metric_array(:,prec_index_num);   
        sharp_right_vec = metric_array(:,sharp_right_index_num);    
        spec_filter = spec_vec>=0&spec_vec<=2*i;
        prec_vec = prec_vec(spec_filter);
        spec_vec = spec_vec(spec_filter);
        sharp_right_vec = sharp_right_vec(spec_filter);
        b_points = boundary(spec_vec,sharp_right_vec,0.95);
    
    
        % identify boundary points    
        sharp_right_b_vec = NaN(1,length(spec_index)-1);
        spec_b_vec = NaN(1,length(spec_index)-1);

        for c = 1:length(spec_index)-1
            spec_filter = find(spec_vec<spec_index(c+1)&spec_vec>=spec_index(c)); 
            if any(spec_filter)
                sharp_right_b_vec(c) = nanmax(sharp_right_vec(spec_filter));        
            end
        end
        if ~high_flag
            master_struct_multi_bs(i).spec_boundary = spec_vec(b_points); 
            master_struct_multi_bs(i).sharp_right_boundary = sharp_right_vec(b_points);                   
        else
            master_struct_multi_bs(i).spec_boundary_high = spec_vec(b_points); 
            master_struct_multi_bs(i).sharp_right_boundary_high = sharp_right_vec(b_points);                   
        end
    end
end
  
%%%%%%%%%%%%%%%%%%%%%%%%5
% Load multi g results
% get list of sweep results files with only 1 genera TF reaction
multi_g_sweep_files = dir([DataPath 'sweep_results_s01_ns00_g0*_cw1_neq.mat']);
multi_g_info_files = dir([DataPath 'sweep_info_s01_ns00_g0*_cw1_neq.mat']);

multi_g_sweep_files_high = dir([DataPath 'sweep_results_s01_ns00_g0*_cw1_neq_cval1000.mat']);
multi_g_info_files_high = dir([DataPath 'sweep_info_s01_ns00_g0*_cw1_neq_cval1000.mat']);

% load
master_struct_multi_g = struct;
for f = 1:length(multi_g_sweep_files)
  
    load([DataPath multi_g_sweep_files(f).name])
    load([DataPath multi_g_info_files(f).name])
    
    master_struct_multi_g(f).sweep_results = sim_results;
    master_struct_multi_g(f).sweep_info = sim_info;
    
    % load high CW results
    load([DataPath multi_g_sweep_files_high(f).name])
    load([DataPath multi_g_info_files_high(f).name])
    
    master_struct_multi_g(f).sweep_results_high = sim_results;
    master_struct_multi_g(f).sweep_info_high = sim_info;
    
end    

rng(321);

for i = 1:length(master_struct_multi_g)
    for high_flag = 0:1
        if ~high_flag
            metric_array = master_struct_multi_g(i).sweep_results.metric_array;          
        else
            metric_array = master_struct_multi_g(i).sweep_results_high.metric_array;          
        end

        % extract vectors        
        spec_vec = metric_array(:,spec_index_num);    
        sharp_right_vec = metric_array(:,sharp_right_index_num);    
        spec_filter = spec_vec>=0&spec_vec<=2*i;
        spec_vec = spec_vec(spec_filter);
        sharp_right_vec = sharp_right_vec(spec_filter);
        b_points = boundary(spec_vec,sharp_right_vec,0.95);


        % identify boundary points    
        sharp_right_b_vec = NaN(1,length(spec_index)-1);
        spec_b_vec = NaN(1,length(spec_index)-1);

        for c = 1:length(spec_index)-1
            spec_filter = find(spec_vec<spec_index(c+1)&spec_vec>=spec_index(c)); 
            if any(spec_filter)
                sharp_right_b_vec(c) = nanmax(sharp_right_vec(spec_filter));        
            end
        end
        
        if ~high_flag
            master_struct_multi_g(i).spec_boundary = spec_vec(b_points); 
            master_struct_multi_g(i).sharp_right_boundary = sharp_right_vec(b_points);  
        else
            master_struct_multi_g(i).spec_boundary_high = spec_vec(b_points); 
            master_struct_multi_g(i).sharp_right_boundary_high = sharp_right_vec(b_points);  
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%5
% Load 2BS multi g results

% get list of sweep results files with only 1 genera TF reaction
multi_g2_sweep_files = dir([DataPath 'sweep_results_s02_ns00_g04_cw1_neq.mat']);
multi_g2_info_files = dir([DataPath 'sweep_info_s02_ns00_g04_cw1_neq.mat']);

multi_g2_sweep_files_high = dir([DataPath 'sweep_results_s02_ns00_g0*_cw1_neq.mat']);
multi_g2_info_files_high = dir([DataPath 'sweep_info_s02_ns00_g0*_cw1_neq.mat']);

% load
master_struct_multi_g2 = struct;
for f = 1:length(multi_g2_sweep_files)
  
    load([DataPath multi_g2_sweep_files(f).name])
    load([DataPath multi_g2_info_files(f).name])
    
    master_struct_multi_g2(f).sweep_results = sim_results;
    master_struct_multi_g2(f).sweep_info = sim_info;
    
    % load high results
    if true%length(multi_g2_sweep_files_high)>=f
        load([DataPath multi_g2_sweep_files_high(f).name])
        load([DataPath multi_g2_info_files_high(f).name])

        master_struct_multi_g2(f).sweep_results_high = sim_results;
        master_struct_multi_g2(f).sweep_info_high = sim_info;
    end
      
end    

rng(321);

for i = 1:length(master_struct_multi_g2)
    for high_flag = 0:1
      
        if ~high_flag
            metric_array = master_struct_multi_g2(i).sweep_results.metric_array;          
        else
            metric_array = master_struct_multi_g2(i).sweep_results_high.metric_array;          
        end

        % extract vectors        
        spec_vec = metric_array(:,spec_index_num);    
        sharp_right_vec = metric_array(:,sharp_right_index_num);    
        spec_filter = spec_vec>=0;%&spec_vec<=2*i;
        spec_vec = spec_vec(spec_filter);
        sharp_right_vec = sharp_right_vec(spec_filter);
        b_points = boundary(spec_vec,sharp_right_vec,0.95);


        % identify boundary points    
        sharp_right_b_vec = NaN(1,length(spec_index)-1);
        spec_b_vec = NaN(1,length(spec_index)-1);

        for c = 1:length(spec_index)-1
            spec_filter = find(spec_vec<spec_index(c+1)&spec_vec>=spec_index(c)); 
            if any(spec_filter)
                sharp_right_b_vec(c) = nanmax(sharp_right_vec(spec_filter));        
            end
        end

        if ~high_flag
            master_struct_multi_g2(i).spec_boundary = spec_vec(b_points); 
            master_struct_multi_g2(i).sharp_right_boundary = sharp_right_vec(b_points);  
        else
            master_struct_multi_g2(i).spec_boundary_high = spec_vec(b_points); 
            master_struct_multi_g2(i).sharp_right_boundary_high = sharp_right_vec(b_points);  
        end
    end
end

%% Plot S0 vs f
close all
markerSize = 75;

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

for high_flag = 0%0:1
    H0_vs_f0_fig = figure;
    hold on

    slb = [];

    % plot 5 BS neq system first
%     if ~high_flag
%         f0_vec_5bs = 10.^master_struct_multi_bs(5).spec_boundary;
%         s0_vec_5bs = master_struct_multi_bs(5).sharp_right_boundary/5;
%     else
%         f0_vec_5bs = 10.^master_struct_multi_bs(5).spec_boundary_high;
%         s0_vec_5bs = master_struct_multi_bs(5).sharp_right_boundary_high/5;
%     end
%     slb(1) = fill(f0_vec_5bs,s0_vec_5bs, cmap_cmb(5,:), 'EdgeColor','k','FaceAlpha',1);

    % 1 BS and multiple conformations     
    for i = fliplr(1:length(master_struct_multi_g)-1)
        if ~high_flag
            f0_vec = 10.^master_struct_multi_g(i).spec_boundary;
            s0_vec = master_struct_multi_g(i).sharp_right_boundary;
        else
            f0_vec = 10.^master_struct_multi_g(i).spec_boundary_high;
            s0_vec = master_struct_multi_g(i).sharp_right_boundary_high;
        end

        slb(end+1) = fill(f0_vec,s0_vec, cmap_pu(2*i,:), 'EdgeColor',brighten(cmap_pu(2*i,:),-0.3),'FaceAlpha',0.5,'LineWidth',3);
    end

    if ~high_flag
        f0_vec = 10.^master_struct_multi_g2(1).spec_boundary;
        s0_vec = master_struct_multi_g2(1).sharp_right_boundary/2;
    else
        f0_vec = 10.^master_struct_multi_g2(1).spec_boundary_high;
        s0_vec = master_struct_multi_g2(1).sharp_right_boundary_high/2;
    end
%     slb(end+1) = fill(f0_vec, s0_vec, cmap_gre(4,:), 'EdgeColor','k','FaceAlpha',1);    

%     slb = [slb(2:end) slb(1)];                
%     legend(slb, 'N_{LC}=2','N_{LC}=3','N_{LC}=4','N_{LC}=5','N_{LC}=4 (N_a=2)','N_{LC}=2 (N_a=5)')
   
    legend(fliplr(slb), 'N_{LC}=2','N_{LC}=3','N_{LC}=4','N_{LC}=5')
    xlabel('specificity gain (f/f_{eq})');
    ylabel('intrinsic sharpness gain (H_0/H_0^{eq})')
    % grid on
    % box on
    set(gca,'FontSize',14)
    set(gca,'xscale','log')
    set(gca,'xtick',[1  10^2 10^4 10^6 10^8])
    set(gca,'xticklabels',{'\alpha^0','\alpha','\alpha^2','\alpha^3','\alpha^4'})
    xlim([1 1e8])
    % ylim([1e-6 1])

    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';
    
    H0_vs_f0_fig.InvertHardcopy = 'off';
    set(gcf,'color','w');
    set(gca,'FontSize',14)
    if ~high_flag
        grid on
        set(gca,'Color',[228,221,209]/255) 
    end
    % ylim([0.25 1.75])

    suffix = '';
    if high_flag
        suffix = '_cw1000';
    end
    saveas(H0_vs_f0_fig,[FigPath 'H0_vs_f0' suffix '.png'])
    saveas(H0_vs_f0_fig,[FigPath 'H0_vs_f0' suffix '.pdf']) 
end
