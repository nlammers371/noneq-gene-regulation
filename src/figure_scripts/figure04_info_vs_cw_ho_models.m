% Plot results for sharpness vs precision for higher order models
clear 
close all
addpath(genpath('../utilities/'))

% %%%%%%%%%%%%%%%%  set relative read and write paths %%%%%%%%%%%%%%%%%%%%
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\';
DataPath = [DropboxFolder  'SweepOutput\sweeps03_info_vs_cw' filesep ];
FigPath = [DropboxFolder '\manuscript\info_vs_cw' filesep];
mkdir(FigPath);

% load dataset with decision time ranges 
load('decision_limit_info.mat','decision_limit_info')

% %%%%%%%%%%%%%%%%  Set plot parameters and constants %%%%%%%%%%%%%%%%%%%%
bit_factor = log2(exp(1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%% Compare "equilibrium" strategy of adding binding sites with 
%%% "Non-equilibrium" strategy of adding locus conformations

% generate decision time range vectors
dt_cw_vec = geomean(vertcat(decision_limit_info.cw_ub,decision_limit_info.cw_lb));
dt_ub_vec = decision_limit_info.n_cycles_ub;
dt_lb_vec = 1e-10*ones(size(dt_cw_vec)); % not worried about lower bound, so exclude from plot
dt_mean_vec = dt_lb_vec;%mean(vertcat(dt_ub_vec,dt_lb_vec));

% combine mous and human
% dt_cw_vec(end-1) = mean(dt_cw_vec(end-1:end));
% dt_cw_vec = dt_cw_vec(1:end-1);
% dt_ub_vec(end-1) = max(dt_ub_vec(end-1:end));
% dt_ub_vec = dt_ub_vec(1:end-1);
% dt_lb_vec(end-1) = min(dt_lb_vec(end-1:end));
% dt_lb_vec = dt_lb_vec(1:end-1);
% 
% dt_mean_vec(end-1) = mean(dt_mean_vec(end-1:end));
% dt_mean_vec = dt_mean_vec(1:end-1);

% exclude Arabadopsis for now
plot_vec = [1 1 0 0 1]==1;

% get metric names for numeric sweeps
[~,~,metric_names_num] = calculateMetricsNumeric_v3([]);

ir_index_num = find(strcmp(metric_names_num,'IR'));
cw_index_num = find(strcmp(metric_names_num,'CW'));
rate_index_num = find(strcmp(metric_names_num,'ProductionRate'));
phi_index_num = find(strcmp(metric_names_num,'Phi'));
tau_index_num = find(strcmp(metric_names_num,'TauCycle'));
inv_dt_index_num = find(strcmp(metric_names_num,'InverseDecisionTime'));
sharp_right_index_num = find(strcmp(metric_names_num,'SharpnessRight'));
spec_index_num = find(strcmp(metric_names_num,'Specificity'));
prec_index_num = find(strcmp(metric_names_num,'Precision'));

cw_index = linspace(-1,6,121);
cw_index = sort([cw_index log10(dt_cw_vec(end))]);

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
    % extract vectors
    metric_array = master_struct_multi_bs(i).sweep_results.metric_array;   
    cw_vec = metric_array(:,cw_index_num);
    ir_vec = metric_array(:,ir_index_num)*bit_factor;            
    dt_vec = 1./metric_array(:,inv_dt_index_num);            
    % identify boundary points
    ir_b_vec = NaN(1,length(cw_index)-1);
    dt_b_vec = NaN(1,length(cw_index)-1);
    for c = 1:length(cw_index)-1
        cw_filter = find(cw_vec<cw_index(c+1)&cw_vec>=cw_index(c));
        [ir_b_vec(c), mi] = nanmax(ir_vec(cw_filter));
        dt_b_vec(c) = dt_vec(cw_filter(mi));
    end
    
    % store    
    master_struct_multi_bs(i).cw_boundary = 10.^cw_index(2:end);
    master_struct_multi_bs(i).ir_boundary = imgaussfilt(ir_b_vec,2);   
    master_struct_multi_bs(i).dt_boundary = imgaussfilt(dt_b_vec,2);   
    
    % do the same for eq systems
    metric_array_eq = master_struct_multi_bs(i).sweep_results_eq.metric_array;   
    cw_vec_eq = metric_array_eq(:,cw_index_num);
    ir_vec_eq = metric_array_eq(:,ir_index_num)*bit_factor;            
    dt_vec_eq = 1./metric_array_eq(:,inv_dt_index_num);   
    % identify boundary points
    ir_b_vec_eq = NaN(1,length(cw_index)-1);
    dt_b_vec_eq = NaN(1,length(cw_index)-1);
    for c = 1:length(cw_index)-1
        cw_filter = find(cw_vec_eq<cw_index(c+1)&cw_vec_eq>=cw_index(c));
        [ir_b_vec_eq(c), mi] = nanmax(ir_vec_eq(cw_filter));
        dt_b_vec_eq(c) = dt_vec_eq(cw_filter(mi));
    end
    
    % store    
    master_struct_multi_bs(i).cw_boundary_eq = 10.^cw_index(2:end);
    master_struct_multi_bs(i).ir_boundary_eq = imgaussfilt(ir_b_vec_eq,2);   
    master_struct_multi_bs(i).dt_boundary_eq = imgaussfilt(dt_b_vec_eq,2);   
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
    % identify boundary points
    ir_b_vec = NaN(1,length(cw_index)-1);
    dt_b_vec = NaN(1,length(cw_index)-1);
    for c = 1:length(cw_index)-1
        cw_filter = find(cw_vec<cw_index(c+1)&cw_vec>=cw_index(c));
        [ir_b_vec(c), mi] = nanmax(ir_vec(cw_filter));
        dt_b_vec(c) = dt_vec(cw_filter(mi));
    end
    
    % store    
    master_struct_multi_g(i).cw_boundary = 10.^cw_index(2:end);
    master_struct_multi_g(i).ir_boundary = imgaussfilt(ir_b_vec,2);  
    master_struct_multi_g(i).dt_boundary = imgaussfilt(dt_b_vec,2);  
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
    % identify boundary points
    ir_b_vec = NaN(1,length(cw_index)-1);
    dt_b_vec = NaN(1,length(cw_index)-1);
    for c = 1:length(cw_index)-1
        cw_filter = find(cw_vec<cw_index(c+1)&cw_vec>=cw_index(c));
        [ir_b_vec(c), mi] = nanmax(ir_vec(cw_filter));
        dt_b_vec(c) = dt_vec(cw_filter(mi));
    end
    
    % store    
    master_struct_multi_g2(i).cw_boundary = 10.^cw_index(2:end);
    master_struct_multi_g2(i).ir_boundary = imgaussfilt(ir_b_vec,2);  
    master_struct_multi_g2(i).dt_boundary = imgaussfilt(dt_b_vec,2);  
end

%% Plot range of achievable information rates
close all

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

       
bottomline = repelem(1e-8,length(cw_index)-1);

% make figure
close all
ir_space_fig_bs = figure('Position',[100 100 560 420]);
hold on
feq = [];
for i = length(master_struct_multi_bs):-1:1
    cw_vec = master_struct_multi_bs(i).cw_boundary_eq;
    ir_vec = master_struct_multi_bs(i).ir_boundary_eq;
    
    fa = 0.75;
    ea = 0.75;
    if i == 1
        fa = 1;
    end
    if i == 1||i==5
        ea = 0.75;
    end
    cm = cmap_cmb(i,:);
   
    feq(end+1) = fill([cw_vec fliplr(cw_vec)],[ir_vec fliplr(bottomline)],...
                          cm,'FaceAlpha',fa,'EdgeAlpha',...
                          ea,'EdgeColor',brighten(cm,-0.5),'LineWidth',3); 
end

% plot Ir noneq bound
% plot(master_struct_multi_bs(end).cw_boundary,master_struct_multi_bs(end).ir_boundary,':','Color',brighten(cmap_cmb(end,:),-0.5),'LineWidth',3)

% legend(feq,'N_a=1','N_a=2','N_a=3','N_a=4','N_a=5')
xlabel('relative wrong factor concentration (w/c)');
ylabel('info. rate (bits per burst cycle)')
grid on
box on
set(gca,'FontSize',18)
% set(gca,'Color',[228,221,209]/255) 
set(gca,'yscale','log')
set(gca,'xscale','log')
set(gca,'xtick',[1  10^2  10^4 ])
xlim([1 1e5])
ylim([1e-6 1])

ax = gca;
% ax.LineWidth = 2;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

ir_space_fig_bs.InvertHardcopy = 'off';
set(gcf,'color','w');
% ylim([0.25 1.75])

saveas(ir_space_fig_bs,[FigPath 'ir_vec_cw_bs.png'])
saveas(ir_space_fig_bs,[FigPath 'ir_vec_cw_bs.pdf']) 

%% Same thing for decision times
topline = repelem(1e10,length(cw_index)-1);
close all
% make figure        
dt_space_fig_bs = figure('Position',[100 100 560 420]);
hold on
feq = [];
for i = length(master_struct_multi_bs):-1:1
    cw_vec = master_struct_multi_bs(i).cw_boundary_eq;
    ir_vec = master_struct_multi_bs(i).dt_boundary_eq;
    
    fa = 0.75;
    ea = 0.75;
    if i == 1
        fa = 1;
    end
    if i == 1||i==5
        ea = 0.75;
    end
    cm = cmap_cmb(i,:);
   
    feq(end+1) = fill([cw_vec fliplr(cw_vec)],[ir_vec fliplr(topline)],...
                          cm,'FaceAlpha',fa,'EdgeAlpha',...
                          ea,'EdgeColor',brighten(cm,-0.5),'LineWidth',3); 
end

% plot IR noneq bound
% plot(master_struct_multi_bs(end).cw_boundary,master_struct_multi_bs(end).dt_boundary,':','Color',brighten(cmap_cmb(end,:),-0.25),'LineWidth',3)

% plot feasible decision time ranges
errorbar(dt_cw_vec(plot_vec),dt_mean_vec(plot_vec),dt_mean_vec(plot_vec)-dt_lb_vec(plot_vec),dt_ub_vec(plot_vec)-dt_mean_vec(plot_vec),...
          '.','LineWidth',1.5','Color','k','CapSize',15)

legend(fliplr(feq),'N_b=1','N_b=2','N_b=3','N_b=4','N_b=5','Location','southeast')
xlabel('relative wrong factor concentration (w/c)');
ylabel('decision time (burst cycles)')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
set(gca,'yscale','log')
set(gca,'xscale','log')
set(gca,'xtick',[1 10 10^2 10^3 10^4 10^5])
xlim([1 1e5])
ylim([1 1e5])

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

dt_space_fig_bs.InvertHardcopy = 'off';
set(gcf,'color','w');
% ylim([0.25 1.75])

saveas(dt_space_fig_bs,[FigPath 'dt_vec_cw_bs.png'])
saveas(dt_space_fig_bs,[FigPath 'dt_vec_cw_bs.pdf']) 

%% Plot range of achievable information rates
close all
ir_space_fig_g = figure('Position',[100 100 560 420]);
hold on

% first plot 2-bs 4-LC system that achieves feasible range
cw_vec = master_struct_multi_g2(3).cw_boundary;
ir_vec = master_struct_multi_g2(3).ir_boundary;
fneq = [];
f2 = fill([cw_vec fliplr(cw_vec)],[ir_vec fliplr(bottomline)],...
                          cmap_cmb(2,:),'FaceAlpha',0.75,'EdgeAlpha',...
                          0.75,'EdgeColor',brighten(cmap_cmb(2,:),-0.5),'LineWidth',3); 
                        
for i = length(master_struct_multi_g)-1:-1:1
    cw_vec = master_struct_multi_g(i).cw_boundary;
    ir_vec = master_struct_multi_g(i).ir_boundary;
    
    fa = 0.75;
    ea = 0.75;
    if i == 1
        fa = 1;
    end
    if i == 1||i==5
        ea = 0.75;
    end
    cm = cmap_pu(2*i,:);
   
    fneq(end+1) = fill([cw_vec fliplr(cw_vec)],[ir_vec fliplr(bottomline)],...
                          cm,'FaceAlpha',fa,'EdgeAlpha',...
                          ea,'EdgeColor',brighten(cm,-0.5),'LineWidth',3); 
end

% plot Ir noneq bound
% plot(master_struct_multi_bs(end).cw_boundary,master_struct_multi_bs(end).ir_boundary,':','Color',brighten(cmap_cmb(end,:),-0.5),'LineWidth',3)

% legend([fliplr(fneq) f2],'N_g=2 (N_a=1)','N_g=3 (N_a=1)','N_g=4 (N_a=1)','N_g=5 (N_a=1)','N_g=4 (N_a=2)','Location','southeast')
xlabel('relative wrong factor concentration (w/c)');
ylabel('info. rate (bits per burst cycle)')

grid on
box on
set(gca,'FontSize',18)
% set(gca,'Color',[228,221,209]/255) 
set(gca,'yscale','log')
set(gca,'xscale','log')
set(gca,'xtick',[1  10^2  10^4])
xlim([1 1e5])
ylim([1e-6 1])

ax = gca;
% ax.LineWidth = 2;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

ir_space_fig_g.InvertHardcopy = 'off';
set(gcf,'color','w');
% ylim([0.25 1.75])

saveas(ir_space_fig_g,[FigPath 'ir_vec_cw_g.png'])
saveas(ir_space_fig_g,[FigPath 'ir_vec_cw_g.pdf']) 

%% Same thing for decision times

topline = repelem(1e10,length(cw_index)-1);
close all

% make figure        
dt_space_fig_g = figure('Position',[100 100 560 420]);
hold on

% first plot 2-bs 4-LC system that achieves feasible range
cw_vec = master_struct_multi_g2(3).cw_boundary;
dt_vec = master_struct_multi_g2(3).dt_boundary;
fneq = [];
f2 = fill([cw_vec fliplr(cw_vec)],[dt_vec fliplr(topline)],...
                          cmap_cmb(2,:),'FaceAlpha',0.75,'EdgeAlpha',...
                          0.75,'EdgeColor',brighten(cmap_cmb(2,:),-0.5),'LineWidth',3); 
                        
for i = length(master_struct_multi_g)-1:-1:1
    cw_vec = master_struct_multi_g(i).cw_boundary;
    ir_vec = master_struct_multi_g(i).dt_boundary;
    
    fa = 0.75;
    ea = 0.75;
    if i == 1
        fa = 1;
        ea = 1;        
    end
    if i==5
        ea = 0.75;
    end
    cm = cmap_pu(2*i,:);
    if false
        cm = cmap2(2,:);
    end
    fneq(end+1) = fill([cw_vec fliplr(cw_vec)],[ir_vec fliplr(topline)],...
                          cm,'FaceAlpha',fa,'EdgeAlpha',...
                          ea,'EdgeColor',brighten(cm,-0.5),'LineWidth',3); 
end

% plot 1bs system at eq
cw_vec = master_struct_multi_bs(1).cw_boundary_eq;
ir_vec = master_struct_multi_bs(1).dt_boundary_eq;

% fill([cw_vec fliplr(cw_vec)],[ir_vec fliplr(topline)],...
%                           cmap2(3,:),'FaceAlpha',0,'EdgeAlpha',...
%                           1,'EdgeColor','k'...brighten(cmap2(3,:),-0.25)
%                           ,'LineWidth',3,'LineStyle',':'); 


% plot IR noneq bound
% plot(master_struct_multi_bs(end).cw_boundary_eq,master_struct_multi_bs(end).dt_boundary_eq,'-','Color',brighten(cmap_cmb(end,:),-0.25),'LineWidth',3)
% plot(master_struct_multi_bs(end).cw_boundary,master_struct_multi_bs(end).dt_boundary,':','Color',brighten(cmap_cmb(end,:),-0.25),'LineWidth',3)

% plot feasible decision time ranges
errorbar(dt_cw_vec(plot_vec),dt_mean_vec(plot_vec),dt_mean_vec(plot_vec)-dt_lb_vec(plot_vec),dt_ub_vec(plot_vec)-dt_mean_vec(plot_vec),...
        '.','LineWidth',1.5','Color','k','CapSize',15)
legend([fliplr(fneq) f2],'N_a=1 (N_b=1)','N_a=2 (N_b=1)','N_a=3 (N_b=1)','N_a=4 (N_b=1)','N_a=3 (N_b=2)','Location','northwest')
xlabel('relative wrong factor concentration (w/c)');
ylabel('decision time (burst cycles)')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
set(gca,'yscale','log')
set(gca,'xscale','log')
set(gca,'xtick',[1 10 10^2 10^3 10^4 10^5])
set(gca,'ytick',[1 10 10^2 10^3 10^4 10^5])
xlim([1 1e5])
ylim([1 1e5])

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

dt_space_fig_g.InvertHardcopy = 'off';
set(gcf,'color','w');
% ylim([0.25 1.75])

saveas(dt_space_fig_g,[FigPath 'dt_vec_cw_g.png'])
saveas(dt_space_fig_g,[FigPath 'dt_vec_cw_g.pdf'])

%% look at trend with BS in equilibrium
close all

[~,cw_i] = min(abs(master_struct_multi_bs(i).cw_boundary_eq-dt_cw_vec(end)));
bs_vec = 1:5;
ir_vec_bs = NaN(1,5);
dt_vec_bs = NaN(1,5);
for i = 1:length(master_struct_multi_bs)
    ir_vec_bs(i) = master_struct_multi_bs(i).ir_boundary_eq(cw_i);
    dt_vec_bs(i) = master_struct_multi_bs(i).dt_boundary_eq(cw_i);
end
bs_vec = [0 bs_vec];
ir_vec_bs = [0 ir_vec_bs];
% dt_vec_bs = [dt_vec_bs];
mdl_ir = fitlm(bs_vec'.^2, ir_vec_bs,'Intercept',false);
mdl_dt = fitlm(bs_vec(2:end)'.^-2, dt_vec_bs,'Intercept',false);

% extrapolate
bs_vec_pd = 0:100;
bs_vec_pd_long = linspace(0.001,100);
ir_vec_eq_pd = predict(mdl_ir,bs_vec_pd'.^2);
dt_vec_eq_pd = predict(mdl_dt,bs_vec_pd'.^-2);
dt_vec_eq_pd_long = predict(mdl_dt,bs_vec_pd_long'.^-2);
% calculate dt coeficient
dtc = log((1-0.32)/.32)*(1-2*0.32);
% dt_vec_eq_pd = dtc./ir_vec_eq_pd;

figure('Position',[100 100 560 420]);
hold on
scatter(bs_vec,ir_vec_bs);
plot(bs_vec_pd,ir_vec_eq_pd)
set(gca,'yscale','log')

ind = find(dt_vec_eq_pd<=dt_ub_vec(end),1);
bs_vec_pd(ind);


% decision time
bs_fig = figure('Position',[100 100 560 420]);
hold on
s = scatter(bs_vec(2:end),dt_vec_bs);
p = plot(bs_vec_pd_long,dt_vec_eq_pd_long');

set(gca,'yscale','log')
ylim([1 1e5])
yl = ylim;

plot(bs_vec_pd,repelem(dt_ub_vec(end),length(dt_vec_eq_pd)),'--','Color','k','LineWidth',2)
plot([bs_vec_pd(ind) bs_vec_pd(ind)],[yl(1) yl(2)],'--','Color','k','LineWidth',2)
ax = gca;

xlabel('number of binding sites')
ylabel('decision time (burst cycles)')
set(gca,'FontSize',14)

legend([s p],'numerical predictions','quadratic extrapolation')
xlim([0 50])
grid on
saveas(bs_fig,[FigPath 'dt_eq_extrap.png'])
saveas(bs_fig,[FigPath 'dt_eq_extrap.pdf'])