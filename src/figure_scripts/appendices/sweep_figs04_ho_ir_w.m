% Plot results for IR vs energy for higher order models
clear 
close all
addpath(genpath('../utilities/'))

% %%%%%%%%%%%%%%%%  set relative read and write paths %%%%%%%%%%%%%%%%%%%%
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\';
DataPath = [DropboxFolder  'SweepOutput\sweeps03B_ir_vs_rate_cw' filesep ];
DataPath_f5 = [DropboxFolder  'SweepOutput\sweeps03_info_vs_cw_cb' filesep ];
FigPath = [DropboxFolder '\manuscript\appendices' filesep 'sweep_algorithm' filesep];
mkdir(FigPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% load numerical sweep results sets MULTI BS
[~,~,metric_names_num] = calculateMetricsNumeric_v3([]);
n_ir = 100;
ir_index = find(strcmp(metric_names_num,'IR'));
w_index = find(strcmp(metric_names_num,'CW'));

%%%%%%%%%
% Multi BS
%%%%%%%%%

% get list of sweep results files with only 1 general TF reaction
ref_sweep_files = dir([DataPath 'sweep_results*']);
ref_info_files = dir([DataPath 'sweep_info*']);

ref_struct = struct; 
% load
for f = 1:length(ref_info_files)
  
    % load files
    load([DataPath  ref_info_files(f).name])
    load([DataPath  ref_sweep_files(f).name])
    
    % add main structure
    ref_struct(f).sweep_results = sim_results;
    ref_struct(f).sweep_info = sim_info;        
    
    % add key info
    ref_struct(f).eq_flag = sim_info.equilibrium_flag;
    ref_struct(f).w = sim_info.cw;
    
    s_ind = strfind(ref_sweep_files(f).name,'_s');
    ref_struct(f).nb = str2double(ref_sweep_files(f).name(s_ind+2:s_ind+3));
    
    g_ind = strfind(ref_sweep_files(f).name,'_g');
    ref_struct(f).nlc = str2double(ref_sweep_files(f).name(g_ind+2:g_ind+3));
    
    % calculate max info
    ref_struct(f).max_ir = nanmax(sim_results.metric_array(:,ir_index))*log2(exp(1));
    
    clear sim_info
    clear sim_results
end            

% generate helper vectors
w_val_vec = [ref_struct.w];
nb_val_vec = [ref_struct.nb];
nlc_val_vec = [ref_struct.nlc];
eq_flag_vec = [ref_struct.eq_flag];
max_ir_vec = [ref_struct.max_ir];

%% Load IR vs. W results used in main text Figure 5
w_index_vec = unique(w_val_vec);

multi_bs_sweep_files_f5 = [dir([DataPath_f5 'sweep_results*g01_cw1_eq*']) ; dir([DataPath_f5 'sweep_results*s05*g01_cw1_neq*'])];
multi_bs_info_files_f5 = [dir([DataPath_f5 'sweep_info*g01_cw1_eq*']) ; dir([DataPath_f5 'sweep_info*s05*g01_cw1_neq*'])];

multi_lc_sweep_files_f5 = dir([DataPath_f5 'sweep_results*s01*_cw1_neq*']);
multi_lc_info_files_f5 = dir([DataPath_f5 'sweep_info*s01*_cw1_neq*']);

master_struct_multi_lc = struct;
master_struct_multi_bs = struct;

for f = 1:length(multi_bs_sweep_files_f5)
  
    % load files
    load([DataPath_f5  multi_bs_info_files_f5(f).name])
    load([DataPath_f5  multi_bs_sweep_files_f5(f).name])
    
    master_struct_multi_bs(f).sweep_results_f5 = sim_results;
    master_struct_multi_bs(f).sweep_info_f5 = sim_info;        
    
    clear sim_results
    clear sim_info
    
    % calculate max IR values for each relevant W level
    w_vec = 10.^master_struct_multi_bs(f).sweep_results_f5.metric_array(:,w_index);
    ir_vec = master_struct_multi_bs(f).sweep_results_f5.metric_array(:,ir_index)*log2(exp(1));
    
    master_struct_multi_bs(f).nb_vec = repelem(min([5,f]),length(w_index_vec));
    master_struct_multi_bs(f).nlc_vec = repelem(1,length(w_index_vec));
    master_struct_multi_bs(f).eq_flag_vec = repelem(master_struct_multi_bs(f).sweep_info_f5.equilibrium_flag,length(w_index_vec));
    ir_max_vec = NaN(size(w_index_vec));
    master_struct_multi_bs(f).w_ref = w_index_vec;
    for i = 1:length(w_index_vec)
        w_filter = w_vec>=w_index_vec(i);
        ir_max_vec(i) = nanmax(ir_vec(w_filter));
    end
    master_struct_multi_bs(f).ir_max_vec = ir_max_vec;
end            

%%%%%%%%%
% Multi-LC
%%%%%%%%%

for f = 1:4
  
    load([DataPath_f5 multi_lc_info_files_f5(f).name])
    load([DataPath_f5  multi_lc_sweep_files_f5(f).name])
    
    master_struct_multi_lc(f).sweep_results_f5 = sim_results;
    master_struct_multi_lc(f).sweep_info_f5 = sim_info;
    
    clear sim_results
    clear sim_info
    
    % calculate max IR values for each relevant W level
    w_vec = 10.^master_struct_multi_lc(f).sweep_results_f5.metric_array(:,w_index);
    ir_vec = master_struct_multi_lc(f).sweep_results_f5.metric_array(:,ir_index)*log2(exp(1));
    
    master_struct_multi_lc(f).nb_vec = repelem(1,length(w_index_vec));
    master_struct_multi_lc(f).nlc_vec = repelem(f,length(w_index_vec));
    master_struct_multi_lc(f).eq_flag_vec = repelem(master_struct_multi_lc(f).sweep_info_f5.equilibrium_flag,length(w_index_vec));
    ir_max_vec = NaN(size(w_index_vec));
    master_struct_multi_lc(f).w_ref = w_index_vec;
    for i = 1:length(w_index_vec)
        w_filter = w_vec>=w_index_vec(i);
        ir_max_vec(i) = nanmax(ir_vec(w_filter));
    end
    master_struct_multi_lc(f).ir_max_vec = ir_max_vec;
end 

% Identify matches and compare IR values

% generate helper vecs
max_ir_vec_full = [[master_struct_multi_lc.ir_max_vec] [master_struct_multi_bs.ir_max_vec]];
w_vec_full = [[master_struct_multi_lc.w_ref] [master_struct_multi_bs.w_ref]];
eq_flag_vec_full = [[master_struct_multi_lc.eq_flag_vec] [master_struct_multi_bs.eq_flag_vec]];
nb_vec_full = [[master_struct_multi_lc.nb_vec] [master_struct_multi_bs.nb_vec]];
nlc_vec_full = [[master_struct_multi_lc.nlc_vec] [master_struct_multi_bs.nlc_vec]];


full_sweep_ir_vals = NaN(size(max_ir_vec));
for i = 1:length(max_ir_vec)
    match_filter = w_vec_full==w_val_vec(i) & eq_flag_vec_full==eq_flag_vec(i) &...
                   nb_vec_full==nb_val_vec(i) & nlc_vec_full==nlc_val_vec(i);
                 
    if sum(match_filter) > 1
        error('too many matches')
    elseif sum(match_filter) == 1
        full_sweep_ir_vals(i) = max_ir_vec_full(match_filter);
    end
end    

%% Make figure
sym_cell = {'o','^','s'};
w_vec = [10 100 1000];
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
        
ir_match_lc = figure;
neq_inds = find(~eq_flag_vec&nb_val_vec==1);

hold on
% plot(linspace(0,0.05),linspace(0,0.05),'-.k','LineWidth',2)
plot(logspace(-5,-1),logspace(-5,-1),'-.k','LineWidth',2)
s = [];
for i = length(neq_inds):-1:1
    nlc = nlc_val_vec(neq_inds(i));
    wi = find(w_val_vec(neq_inds(i))==w_vec);
    s(nlc) = scatter(max_ir_vec(neq_inds(i)),full_sweep_ir_vals(neq_inds(i)),75,sym_cell{wi},...
            'MarkerFaceColor',cmap_pu(2+nlc,:),'MarkerEdgeColor','k','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);
end


xlabel('maximum IR value (IR vs. r)');
ylabel('maximum IR value (IR vs. w)');

grid on
set(gca,'FontSize',14)
legend(s,'N_{A}=1','N_{A}=2','N_{A}=3','N_{A}=4','Location','southeast')
set(gca,'Color',[228,221,209]/255) 

set(gca,'yscale','log')
set(gca,'xscale','log')

% ylim([0 0.05])
% xlim([0 0.05])
ylim([3e-4 0.1])
xlim([3e-4 0.1])

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';


ir_match_lc.InvertHardcopy = 'off';
set(gcf,'color','w');      
saveas(ir_match_lc,[FigPath 'ir_w_scatter_NA.png'])   
saveas(ir_match_lc,[FigPath 'ir_w_scatter_NA.pdf']) 

%% 
        
ir_match_bs = figure;
eq_inds = find(eq_flag_vec|nb_val_vec==5);

hold on
plot(logspace(-5,0),logspace(-5,0),'-.k','LineWidth',2)
% plot(linspace(0,0.25),linspace(0,0.25),'-.k','LineWidth',2)
s = [];
for i = 1:length(eq_inds)%:-1:1
    nb = nb_val_vec(eq_inds(i));
    eq_flag = eq_flag_vec(eq_inds(i));
    wi = find(w_val_vec(eq_inds(i))==w_vec);
    if eq_flag
        sb = scatter(max_ir_vec(eq_inds(i)),full_sweep_ir_vals(eq_inds(i)),75,sym_cell{wi},...
                'MarkerFaceColor',cmap_full(nb,:),'MarkerEdgeColor','k','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);
        if w_val_vec(eq_inds(i))==10
            s(end+1) = sb;
        end
    else
        sb = scatter(max_ir_vec(eq_inds(i)),full_sweep_ir_vals(eq_inds(i)),100,sym_cell{wi},...
                'MarkerEdgeColor',brighten(cmap_full(nb,:),-0.5),'LineWidth',2,'MarkerEdgeAlpha',1);
    end
end

xlabel('maximum IR value (IR vs. r)');
ylabel('maximum IR value (IR vs. w)');

set(gca,'yscale','log')
set(gca,'xscale','log')

grid on
set(gca,'FontSize',14)
legend(s,'N_{B}=1','N_{B}=2','N_{B}=3','N_{B}=4','N_{B}=5','Location','southeast')
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

% ylim([0 0.07])
% xlim([0 0.07])

ir_match_bs.InvertHardcopy = 'off';
set(gcf,'color','w');      
saveas(ir_match_bs,[FigPath 'ir_w_scatter_NB.png'])   
saveas(ir_match_bs,[FigPath 'ir_w_scatter_NB.pdf']) 

%% Make illustrative figure of rate vs IR

nb = 3;
eq = 1;
nlc = 1;

plot_inds = find(nb_val_vec==nb & nlc_val_vec==nlc & eq_flag_vec == eq);

ir_vs_r_fig = figure;
cmap = brewermap(5,'Reds');
hold on
for p = 1:length(plot_inds)
    scatter(ref_struct(plot_inds(p)).sweep_results.metric_array(:,1),ref_struct(plot_inds(p)).sweep_results.metric_array(:,5),75,...
      'MarkerFaceColor',cmap(p+1,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
end    

xlabel('transcription rate (r)');
ylabel('information per burst (bits)');

set(gca,'yscale','log')
ylim([1e-5 0.5])

grid on
set(gca,'FontSize',14)
legend('w/c=10','w/c=10^2','w/c=10^3','Location','northwest')
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

% ylim([0 0.07])
% xlim([0 0.07])

ir_vs_r_fig.InvertHardcopy = 'off';
set(gcf,'color','w');      

saveas(ir_vs_r_fig,[FigPath 'ir_r_example_BS3.png'])  

%% Check convergence of IR vs. r runs

cv_frac_vec = NaN(1,length(ref_struct));
area_frac_vec = NaN(1,length(ref_struct));
area_array = NaN(250,1000,length(ref_struct));

for r = 1:length(ref_struct)
    cv_frac_vec(r) = mean(ref_struct(r).sweep_info.convergence_flag_vec);
    
    area_vec = ref_struct(r).sweep_results.areaVec;
    area_id_vec = ref_struct(r).sweep_results.areaIDVec;
    a_array_temp = area_array(:,:,r);
    a_index = unique(area_id_vec);
    for a = 1:length(a_index)
        a_vec = area_vec(area_id_vec==a_index(a));
        a_array_temp(1:length(a_vec),a) = a_vec;
        a_array_temp(length(a_vec)+1:end,a) = a_vec(end);
    end
    a_array_temp = a_array_temp ./ nanmax(a_array_temp(:));
    area_array(:,:,r) = a_array_temp;
    area_frac_vec(r) = nansum(a_array_temp(end,:)>=0.95)/sum(~isnan(a_array_temp(end,:)));
end    
    
mean(area_frac_vec)    
mean(cv_frac_vec)    