% Plot results for IR vs energy for higher order models
clear 
close all
addpath(genpath('../utilities/'))

% %%%%%%%%%%%%%%%%  set relative read and write paths %%%%%%%%%%%%%%%%%%%%
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\';
DataPath_ir_phi = [DropboxFolder  'SweepOutput\sweeps01_info_vs_energy' filesep ];
DataPath_SP = [DropboxFolder  'SweepOutput\sweeps02_sharpness_vs_precision' filesep];
DataPath_ir_w = [DropboxFolder  'SweepOutput\sweeps03_info_vs_cw_v2' filesep];
FigPath = [DropboxFolder '\manuscript\appendices' filesep 'sweep_algorithm' filesep];
mkdir(FigPath);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% load numerical sweep results sets MULTI BS

%%%%%%%%%
%% S VS P 
%%%%%%%%%

% get list of sweep results files with only 1 genera TF reaction
multi_bs_sweep_files_s_p = dir([DataPath_SP 'sweep_results*g01*neq*']);
multi_bs_info_files_s_p = dir([DataPath_SP 'sweep_info*g01*neq*']);

% load
master_struct_multi_bs = struct;
for f = 1:length(multi_bs_sweep_files_s_p)
  
    load([DataPath_SP multi_bs_sweep_files_s_p(f).name])
    load([DataPath_SP multi_bs_info_files_s_p(f).name])
    
    master_struct_multi_bs(f).sweep_results_sp = sim_results;
    master_struct_multi_bs(f).sweep_info_sp = sim_info;
end    

% add eq files
multi_bs_sweep_files_s_p_eq = dir([DataPath_SP 'sweep_results*g01_cw0_eq*']);
multi_bs_info_files_s_p_eq = dir([DataPath_SP 'sweep_info*g01_cw0_eq*']);
for f = 1:length(multi_bs_sweep_files_s_p_eq)
  
    load([DataPath_SP multi_bs_sweep_files_s_p_eq(f).name])
    load([DataPath_SP multi_bs_info_files_s_p_eq(f).name])
    
    master_struct_multi_bs(f).sweep_results_sp_eq = sim_results;
    master_struct_multi_bs(f).sweep_info_sp_eq = sim_info;
end    

%%%%%%%%%
% IR VS PHI
%%%%%%%%%

% get list of sweep results files with only 1 genera TF reaction
multi_bs_sweep_files_ir = dir([DataPath_ir_phi 'sweep_results*g01*']);
multi_bs_info_files_ir = dir([DataPath_ir_phi 'sweep_info*g01*']);

% load
% master_struct_multi_bs = struct;
for f = 1:length(multi_bs_sweep_files_ir)
  
    load([DataPath_ir_phi  multi_bs_sweep_files_ir(f).name])
    load([DataPath_ir_phi  multi_bs_info_files_ir(f).name])
    
    master_struct_multi_bs(f).sweep_results_ir_phi = sim_results;
    master_struct_multi_bs(f).sweep_info_ir_phi = sim_info;
end            

% add eq files
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load numerical sweep results sets MULTI LC

%%%%%%%%%
% S VS P 
%%%%%%%%%

% get list of sweep results files with only 1 genera TF reaction
multi_lc_sweep_files_s_p = dir([DataPath_SP 'sweep_results*s01*neq*']);
multi_lc_info_files_s_p = dir([DataPath_SP 'sweep_info*s01*neq*']);

% load
master_struct_multi_lc = struct;
for f = 1:length(multi_lc_sweep_files_s_p)-1
  
    load([DataPath_SP multi_lc_sweep_files_s_p(f).name])
    load([DataPath_SP multi_lc_info_files_s_p(f).name])
    
    master_struct_multi_lc(f).sweep_results_sp = sim_results;
    master_struct_multi_lc(f).sweep_info_sp = sim_info;
end    

% add eq files
multi_lc_sweep_files_s_p_eq = dir([DataPath_SP 'sweep_results*g01_cw0_eq*']);
multi_lc_info_files_s_p_eq = dir([DataPath_SP 'sweep_info*g01_cw0_eq*']);
for f = 1:length(multi_lc_sweep_files_s_p_eq)-1
  
    load([DataPath_SP multi_lc_sweep_files_s_p_eq(f).name])
    load([DataPath_SP multi_lc_info_files_s_p_eq(f).name])
    
    master_struct_multi_lc(f).sweep_results_sp_eq = sim_results;
    master_struct_multi_lc(f).sweep_info_sp_eq = sim_info;
end 
%%%%%%%%%
% IR VS PHI
%%%%%%%%%

% get list of sweep results files with only 1 genera TF reaction
multi_lc_sweep_files_ir = dir([DataPath_ir_phi 'sweep_results*s01*']);
multi_lc_info_files_ir = dir([DataPath_ir_phi 'sweep_info*s01*']);

% load
% master_struct_multi_lc = struct;
for f = 1:length(multi_lc_sweep_files_ir)
  
    load([DataPath_ir_phi  multi_lc_sweep_files_ir(f).name])
    load([DataPath_ir_phi  multi_lc_info_files_ir(f).name])
    
    master_struct_multi_lc(f).sweep_results_ir_phi = sim_results;
    master_struct_multi_lc(f).sweep_info_ir_phi = sim_info;
end    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic convergence plots for S vs. P, for which convergence is well-defined

cv_struct_multi_bs = struct;

% S vs. P
for i = 1:length(master_struct_multi_bs)
    % NEQ
    sp_area_id_vec = master_struct_multi_bs(i).sweep_results_sp.areaIDVec;
    id_index = unique(sp_area_id_vec);
    sp_area_vec = master_struct_multi_bs(i).sweep_results_sp.areaVec;
        
    % initialize area array
    area_array = NaN(master_struct_multi_bs(i).sweep_info_sp.n_iters_max,length(id_index));
    area_array_ext = area_array;
    for j = 1:length(id_index)
        a_vec = sp_area_vec(sp_area_id_vec==id_index(j));
        area_array(1:length(a_vec),j) = a_vec;
        area_array_ext(1:length(a_vec),j) = a_vec;
        area_array_ext(length(a_vec)+1:end,j) = a_vec(end);
    end
    cv_struct_multi_bs(i).area_array_sp_neq = area_array;
    area_array_norm = area_array ./ nanmax(area_array(:));
    cv_struct_multi_bs(i).area_array_sp_neq_norm = area_array_norm;
    cv_struct_multi_bs(i).area_array_sp_neq_ext = area_array_ext;
    area_array_ext_norm = area_array_ext ./ nanmax(area_array_ext(:));
    cv_struct_multi_bs(i).area_array_sp_neq_norm_ext = area_array_ext_norm;
    cv_struct_multi_bs(i).cv_frac_sp_neq = mean(master_struct_multi_bs(i).sweep_results_sp.convergence_flags);  
    cv_struct_multi_bs(i).cv_flags_neq = master_struct_multi_bs(i).sweep_results_sp.convergence_flags;
    
    % EQ
    sp_area_id_vec = master_struct_multi_bs(i).sweep_results_sp_eq.areaIDVec;
    id_index = unique(sp_area_id_vec);
    sp_area_vec = master_struct_multi_bs(i).sweep_results_sp_eq.areaVec;
        
    % initialize area array
    area_array = NaN(master_struct_multi_bs(i).sweep_info_sp_eq.n_iters_max,length(id_index));
    area_array_ext = area_array;
    for j = 1:length(id_index)
        a_vec = sp_area_vec(sp_area_id_vec==id_index(j));
        area_array(1:length(a_vec),j) = a_vec;
        area_array_ext(1:length(a_vec),j) = a_vec;
        area_array_ext(length(a_vec)+1:end,j) = a_vec(end);
    end
    cv_struct_multi_bs(i).area_array_sp_eq = area_array;
    area_array_norm = area_array ./ nanmax(area_array(:));
    cv_struct_multi_bs(i).area_array_sp_eq_norm = area_array_norm;
    cv_struct_multi_bs(i).area_array_sp_eq_ext = area_array_ext;
    area_array_ext_norm = area_array_ext ./ nanmax(area_array_ext(:));
    cv_struct_multi_bs(i).area_array_sp_eq_norm_ext = area_array_ext_norm;
    cv_struct_multi_bs(i).cv_frac_sp_eq = mean(master_struct_multi_bs(i).sweep_results_sp_eq.convergence_flags);  
    cv_struct_multi_bs(i).cv_flags_eq = master_struct_multi_bs(i).sweep_results_sp_eq.convergence_flags;
end   
 
% LC S vs. P
cv_struct_multi_lc = struct;
for i = 1:length(master_struct_multi_lc)
    % NEQ
    sp_area_id_vec = master_struct_multi_lc(i).sweep_results_sp.areaIDVec;
    id_index = unique(sp_area_id_vec);
    sp_area_vec = master_struct_multi_lc(i).sweep_results_sp.areaVec;
        
    % initialize area array
    area_array = NaN(master_struct_multi_lc(i).sweep_info_sp.n_iters_max,length(id_index));
    area_array_ext = area_array;
    for j = 1:length(id_index)
        a_vec = sp_area_vec(sp_area_id_vec==id_index(j));
        area_array(1:length(a_vec),j) = a_vec;
        area_array_ext(1:length(a_vec),j) = a_vec;
        area_array_ext(length(a_vec)+1:end,j) = a_vec(end);
    end
    cv_struct_multi_lc(i).area_array_sp_neq = area_array;
    area_array_norm = area_array ./ nanmax(area_array(:));
    cv_struct_multi_lc(i).area_array_sp_neq_norm = area_array_norm;
    cv_struct_multi_lc(i).area_array_sp_neq_ext = area_array_ext;
    area_array_ext_norm = area_array_ext ./ nanmax(area_array_ext(:));
    cv_struct_multi_lc(i).area_array_sp_neq_norm_ext = area_array_ext_norm;
    cv_struct_multi_lc(i).cv_frac_sp_neq = mean(master_struct_multi_lc(i).sweep_results_sp.convergence_flags);  
    cv_struct_multi_lc(i).cv_flags_neq = master_struct_multi_lc(i).sweep_results_sp.convergence_flags;
    
    % EQ
    sp_area_id_vec = master_struct_multi_lc(i).sweep_results_sp_eq.areaIDVec;
    id_index = unique(sp_area_id_vec);
    sp_area_vec = master_struct_multi_lc(i).sweep_results_sp_eq.areaVec;
        
    % initialize area array
    area_array = NaN(master_struct_multi_lc(i).sweep_info_sp_eq.n_iters_max,length(id_index));
    area_array_ext = area_array;
    for j = 1:length(id_index)
        a_vec = sp_area_vec(sp_area_id_vec==id_index(j));
        area_array(1:length(a_vec),j) = a_vec;
        area_array_ext(1:length(a_vec),j) = a_vec;
        area_array_ext(length(a_vec)+1:end,j) = a_vec(end);
    end
    cv_struct_multi_lc(i).area_array_sp_eq = area_array;
    area_array_norm = area_array ./ nanmax(area_array(:));
    cv_struct_multi_lc(i).area_array_sp_eq_norm = area_array_norm;
    cv_struct_multi_lc(i).area_array_sp_eq_ext = area_array_ext;
    area_array_ext_norm = area_array_ext ./ nanmax(area_array_ext(:));
    cv_struct_multi_lc(i).area_array_sp_eq_norm_ext = area_array_ext_norm;
    cv_struct_multi_lc(i).cv_frac_sp_eq = mean(master_struct_multi_lc(i).sweep_results_sp_eq.convergence_flags);  
    cv_struct_multi_lc(i).cv_flags_eq = master_struct_multi_lc(i).sweep_results_sp_eq.convergence_flags;
end

%%
% lc_plot_ind = 4;

lc_sp_fig = figure;
hold on
% normalize array
a_array_norm_neq_5 = cv_struct_multi_lc(4).area_array_sp_neq_norm;
a_array_norm_neq_2 = cv_struct_multi_lc(1).area_array_sp_neq_norm;
n_iters = size(a_array_norm_neq_2,2);
n_steps_max = size(a_array_norm_neq_2,1);

cmap = brewermap(n_iters*6,'Purples');

plot(0:n_steps_max,repelem(1,n_steps_max+1),'-.k','LineWidth',2);

for a = randsample(1:n_iters,n_iters,false)
    p1 = plot(0:n_steps_max,[0 ; a_array_norm_neq_5(:,a)],'Color',[cmap(5*n_iters+a,:) 0.75]);    
end    
for a = randsample(1:n_iters,n_iters,false)
    p2 = plot(0:n_steps_max,[0 ; a_array_norm_neq_2(:,a)],'Color',[cmap(2*n_iters+a,:) 0.75]);
end    

xlabel('sweep iteration');
ylabel('Area')

legend([p2 p1], 'N_{LC}=2','N_{LC}=5','Location','southeast')

grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

ylim([0 1.1])

lc_sp_fig.InvertHardcopy = 'off';
set(gcf,'color','w');      
saveas(lc_sp_fig,[FigPath 'area_plot_lc_sp.png'])   
saveas(lc_sp_fig,[FigPath 'area_plot_lc_sp.pdf'])   

% Define colormaps for use throughout
cmap_pu = brewermap(8,'Purples');
cmap_rd = brewermap(8,'Reds');
cmap_bu = brewermap(8,'Blues');
cmap_gre = brewermap(8,'Greens');
cmap_gra = brewermap(8,'Greys');

% close all
color_ind = 5;
cmap_full = [cmap_pu(3,:); cmap_gre(color_ind,:); cmap_rd(color_ind,:); ...
          cmap_bu(color_ind,:) ; cmap_gra(color_ind,:)];
        
        
bs_sp_fig = figure;
hold on
% normalize array
a_array_norm_neq_5 = cv_struct_multi_bs(5).area_array_sp_neq_norm;
a_array_norm_neq_1 = cv_struct_multi_bs(1).area_array_sp_neq_norm;

cmap_gra = brewermap(6*n_iters,'Greys');

plot(0:n_steps_max,repelem(1,n_steps_max+1),'-.k','LineWidth',2);

for a = randsample(1:n_iters,n_iters,false)
    p1 = plot(0:n_steps_max,[0 ; a_array_norm_neq_5(:,a)],'Color',[cmap_gra(4*n_iters+a,:) 0.25]);    
end    
for a = randsample(1:n_iters,n_iters,false)
    p2 = plot(0:n_steps_max,[0 ; a_array_norm_neq_1(:,a)],'Color',[cmap(2*n_iters+a,:) 1]);
end    

xlabel('sweep iteration');
ylabel('Area')

legend([p2 p1], 'N_{B}=1','N_{B}=5','Location','southeast')

grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

ylim([0 1.1])

bs_sp_fig.InvertHardcopy = 'off';
set(gcf,'color','w');      
saveas(bs_sp_fig,[FigPath 'area_plot_bs_sp.png'])    
saveas(bs_sp_fig,[FigPath 'area_plot_bs_sp.pdf'])  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot fraction of runs converging
close all

lc_cv_frac_eq = [cv_struct_multi_lc.cv_frac_sp_eq];
lc_cv_frac_neq = [cv_struct_multi_lc.cv_frac_sp_neq];

lc_sp_cv = figure;
hold on
cmap_pu = brewermap(8,'Purples');

plot(1:4,lc_cv_frac_neq,'Color','k','LineWidth',2)
for b = 1:4
    scatter(b,lc_cv_frac_neq(b),100,'o','MarkerFaceColor',cmap_pu(2+b,:),'MarkerEdgeColor','k')
end 

plot(1:4,lc_cv_frac_eq,'Color','k','LineWidth',2)
for b = 1:4
    scatter(b,lc_cv_frac_eq(b),75,'s','MarkerFaceColor',cmap_pu(2+b,:),'MarkerEdgeColor','k')
end    

sneq = scatter(-10,-10,100,'o','MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor','k');
seq = scatter(-10,-10,100,'s','MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor','k');

xlabel('number of locus conformations (N_{LC})');
ylabel('fraction of runs converged')


grid on
set(gca,'FontSize',14)
legend([seq sneq], 'equilibrium','non-equilibrium','Location','southwest')
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
ylim([0 1.1])
set(gca,'xtick',1:5)

lc_sp_cv.InvertHardcopy = 'off';
set(gcf,'color','w');      
saveas(lc_sp_cv,[FigPath 'cv_frac_lc_sp.png'])   
saveas(lc_sp_cv,[FigPath 'cv_frac_lc_sp.pdf']) 


bs_cv_frac_eq = [cv_struct_multi_bs.cv_frac_sp_eq];
bs_cv_frac_neq = [cv_struct_multi_bs.cv_frac_sp_neq];

bs_sp_cv = figure;
hold on

plot(1:5,bs_cv_frac_neq,'Color','k','LineWidth',2)
for b = 1:5
    scatter(b,bs_cv_frac_neq(b),100,'o','MarkerFaceColor',cmap_full(b,:),'MarkerEdgeColor','k')
end 

plot(1:5,bs_cv_frac_eq,'Color','k','LineWidth',2)
for b = 1:5
    scatter(b,bs_cv_frac_eq(b),75,'s','MarkerFaceColor',cmap_full(b,:),'MarkerEdgeColor','k')
end    

sneq = scatter(-10,-10,100,'o','MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor','k');
seq = scatter(-10,-10,100,'s','MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor','k');

xlabel('number of binding sites (N_B)');
ylabel('fraction of runs converged')

grid on
set(gca,'FontSize',14)
legend([seq sneq], 'equilibrium','non-equilibrium','Location','southwest')
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

ylim([0 1.1])
set(gca,'xtick',1:5)

bs_sp_cv.InvertHardcopy = 'off';
set(gcf,'color','w');      
saveas(bs_sp_cv,[FigPath 'cv_frac_bs_sp.png'])   
saveas(bs_sp_cv,[FigPath 'cv_frac_bs_sp.pdf'])  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot fraction of runs reaching 95% of max
close all
cv_frac = 0.95;
n_iters = 500;

eq_frac_array_lc = NaN(1,length(cv_struct_multi_lc));
neq_frac_array_lc = NaN(1,length(cv_struct_multi_lc));
for n = 1:length(cv_struct_multi_lc)
    eq_frac_array_lc(n) = nanmean(cv_struct_multi_lc(n).area_array_sp_eq_norm_ext(end,:)>=cv_frac);
    neq_frac_array_lc(n) = nanmean(cv_struct_multi_lc(n).area_array_sp_neq_norm_ext(end,:)>=cv_frac);
end    


eq_frac_array_bs = NaN(1,length(cv_struct_multi_bs));
neq_frac_array_bs = NaN(1,length(cv_struct_multi_bs));
for n = 1:length(cv_struct_multi_bs)
    eq_frac_array_bs(n) = nanmean(cv_struct_multi_bs(n).area_array_sp_eq_norm_ext(end,:)>=cv_frac);
    neq_frac_array_bs(n) = nanmean(cv_struct_multi_bs(n).area_array_sp_neq_norm_ext(end,:)>=cv_frac);
end    

lc_sp_cv = figure;
hold on
cmap_pu = brewermap(8,'Purples');

plot(1:4,neq_frac_array_lc,'Color','k','LineWidth',2)
for b = 1:4
    scatter(b,neq_frac_array_lc(b),100,'o','MarkerFaceColor',cmap_pu(2+b,:),'MarkerEdgeColor','k')
end 

plot(1:4,eq_frac_array_lc,'Color','k','LineWidth',2)
for b = 1:4
    scatter(b,eq_frac_array_lc(b),75,'s','MarkerFaceColor',cmap_pu(2+b,:),'MarkerEdgeColor','k')
end    

sneq = scatter(-10,-10,100,'o','MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor','k');
seq = scatter(-10,-10,100,'s','MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor','k');

xlabel('number of locus conformations (N_{LC})');
ylabel('fraction of runs within 95% of max area')

grid on
set(gca,'FontSize',14)
legend([seq sneq], 'equilibrium','non-equilibrium','Location','southwest')
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

ylim([0 1.1])
set(gca,'xtick',1:5)

lc_sp_cv.InvertHardcopy = 'off';
set(gcf,'color','w');      
saveas(lc_sp_cv,[FigPath 'area_cv_frac_lc_sp.png'])   
saveas(lc_sp_cv,[FigPath 'area_cv_frac_lc_sp.pdf']) 


bs_sp_cv = figure;
hold on

plot(1:5,neq_frac_array_bs,'Color','k','LineWidth',2)
for b = 1:5
    scatter(b,neq_frac_array_bs(b),100,'o','MarkerFaceColor',cmap_full(b,:),'MarkerEdgeColor','k')
end 

plot(1:5,eq_frac_array_bs,'Color','k','LineWidth',2)
for b = 1:5
    scatter(b,eq_frac_array_bs(b),75,'s','MarkerFaceColor',cmap_full(b,:),'MarkerEdgeColor','k')
end    

sneq = scatter(-10,-10,100,'o','MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor','k');
seq = scatter(-10,-10,100,'s','MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor','k');

xlabel('number of binding sites (N_B)');
ylabel('fraction of runs within 95% of max area')

grid on
set(gca,'FontSize',14)
legend([seq sneq], 'equilibrium','non-equilibrium','Location','southwest')
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

ylim([0 1.1])
set(gca,'xtick',1:5)

bs_sp_cv.InvertHardcopy = 'off';
set(gcf,'color','w');      
saveas(bs_sp_cv,[FigPath 'area_cv_frac_bs_sp.png'])   
saveas(bs_sp_cv,[FigPath 'area_cv_frac_bs_sp.pdf'])  


%% Cross-validate max IR values between S vs P and IR vs. Phi sweeps
% because there is no natural boundary along the Phi axis, total area is
% not a good indication of model convergence. Thus, to check that IR vs.
% Phi sweeps are enumerating the full relevant parameter space, we will
% compare maximum IR values for these sweeps with those from S vs P sweeps
% (for which we deomonstrated convergence above). 
[~,~,metric_names_num] = calculateMetricsNumeric_v3([]);
ir_index_num = find(strcmp(metric_names_num,'IR'));

% calculate maximum IR values for all relevant model types 
for i = 1:length(master_struct_multi_lc)
    % NEQ SP
    metric_array_sp = master_struct_multi_lc(i).sweep_results_sp.metric_array;         
    ir_vec = metric_array_sp(:,ir_index_num)*log2(exp(1));
    cv_struct_multi_lc(i).ir_max_sp = nanmax(ir_vec);     
    
    % NEQ IR-PHI
    metric_array_ip = master_struct_multi_lc(i).sweep_results_ir_phi.metric_array;         
    ir_vec = metric_array_ip(:,ir_index_num)*log2(exp(1));
    cv_struct_multi_lc(i).ir_max_ip = nanmax(ir_vec); 
end    

for i = 1:length(master_struct_multi_bs)
    % NEQ SP
    metric_array_sp = master_struct_multi_bs(i).sweep_results_sp.metric_array;         
    ir_vec = metric_array_sp(:,ir_index_num)*log2(exp(1));
    cv_struct_multi_bs(i).ir_max_sp = nanmax(ir_vec); 
    cv_struct_multi_bs(i).cv_frac_sp = mean(master_struct_multi_bs(i).sweep_results_sp.convergence_flags);
    
    % NEQ IR-PHI
    metric_array_ip = master_struct_multi_bs(i).sweep_results_ir_phi.metric_array;         
    ir_vec = metric_array_ip(:,ir_index_num)*log2(exp(1));
    cv_struct_multi_bs(i).ir_max_ip = nanmax(ir_vec); 
    
end    


ir_max_sp_bs = [cv_struct_multi_bs.ir_max_sp];
ir_max_ip_bs = [cv_struct_multi_bs.ir_max_ip];

ir_max_sp_lc = [cv_struct_multi_lc.ir_max_sp];
ir_max_ip_lc = [cv_struct_multi_lc.ir_max_ip];

ir_match = figure;
hold on
cmap_rd = brewermap(3,'Reds');
cmap_pu = brewermap(3,'Purples');

scatter(ir_max_sp_bs,ir_max_ip_bs,100,'MarkerFaceColor',cmap_rd(2,:),'MarkerEdgeColor','k')
scatter(ir_max_sp_lc,ir_max_ip_lc,100,'s','MarkerFaceColor',cmap_pu(2,:),'MarkerEdgeColor','k')

xlabel('maximum IR value (S vs. P sweeps)');
ylabel('maximum IR value (IR vs. \Phi sweeps)');

grid on
set(gca,'FontSize',14)
legend('multi-binding site','multi-locus conformation','Location','southeast')
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';


ir_match.InvertHardcopy = 'off';
set(gcf,'color','w');      
saveas(ir_match,[FigPath 'ir_scatter.png'])   
saveas(ir_match,[FigPath 'ir_scatter.pdf']) 
