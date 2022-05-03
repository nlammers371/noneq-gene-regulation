% Plot results for sharpness vs precision for higher order models
clear 
close all
addpath(genpath('../utilities/'))

% %%%%%%%%%%%%%%%%  set relative read and write paths %%%%%%%%%%%%%%%%%%%%
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\';
DataPath = [DropboxFolder  'SweepOutput\sweeps03_info_vs_cw' filesep ];
FigPath = [DropboxFolder '\manuscript\writeup' filesep];
mkdir(FigPath);

% load dataset with decision time ranges 
load('decision_limit_info.mat','decision_limit_info')

% %%%%%%%%%%%%%%%%  Set plot parameters and constants %%%%%%%%%%%%%%%%%%%%
n_plot = 3e3; % number of points to plot
n_samp = 5e4; % number of points to plot
markerAlpha = 0.5; % marker transparency
bit_factor = log2(exp(1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%% Compare "equilibrium" strategy of adding binding sites with 
%%% "Non-equilibrium" strategy of adding locus conformations

% generate decision time range vectors
dt_cw_vec = geomean(vertcat(decision_limit_info.cw_ub,decision_limit_info.cw_lb));
dt_ub_vec = decision_limit_info.n_cycles_ub;
dt_lb_vec = max(vertcat(decision_limit_info.n_cycles_lb,ones(size(dt_cw_vec)))); % cap lower limit at 1 cycle
dt_mean_vec = mean(vertcat(dt_ub_vec,dt_lb_vec));

% combine mous and human
dt_cw_vec(end-1) = mean(dt_cw_vec(end-1:end));
dt_cw_vec = dt_cw_vec(1:end-1);
dt_ub_vec(end-1) = max(dt_ub_vec(end-1:end));
dt_ub_vec = dt_ub_vec(1:end-1);
dt_lb_vec(end-1) = min(dt_lb_vec(end-1:end));
dt_lb_vec = dt_lb_vec(1:end-1);

dt_mean_vec(end-1) = mean(dt_mean_vec(end-1:end));
dt_mean_vec = dt_mean_vec(1:end-1);

% exclude Arabadopsis for now
plot_vec = [1 1 0 1]==1;

% get metric names for numeric sweeps
[~,~,metric_names_num] = calculateMetricsNumeric_v3([]);

ir_index_num = find(strcmp(metric_names_num,'IR'));
cw_index_num = find(strcmp(metric_names_num,'CW'));
rate_index_num = find(strcmp(metric_names_num,'ProductionRate'));
phi_index_num = find(strcmp(metric_names_num,'Phi'));
tau_index_num = find(strcmp(metric_names_num,'TauCycle'));
inv_dt_index_num = find(strcmp(metric_names_num,'InverseDecisionTime'));


cw_index = linspace(-1,6,121);
cw_index = sort([cw_index log10(dt_cw_vec(end))]);

%%%%%%%%%%%%%%%%%%%%%%%%5
%%% Load multi bs results
% get list of sweep results files with only 1 genera TF reaction
multi_bs_sweep_files_eq = dir([DataPath 'sweep_results*g01*eq*']);
multi_bs_info_files_eq = dir([DataPath 'sweep_info*g01*eq*']);


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
    dt_vec = metric_array(:,ir_index_num)*bit_factor;    
    dt_vec = 1./metric_array(:,inv_dt_index_num); 
    % identify boundary points
    ir_b_vec = NaN(1,length(cw_index)-1);
    dt_b_vec = NaN(1,length(cw_index)-1);
    for c = 1:length(cw_index)-1
        cw_filter = find(cw_vec<cw_index(c+1)&cw_vec>=cw_index(c));
        [ir_b_vec(c), mi] = nanmax(dt_vec(cw_filter));
        dt_b_vec(c) = dt_vec(cw_filter(mi));
    end
    
    % store    
    master_struct_multi_g(i).cw_boundary = 10.^cw_index(2:end);
    master_struct_multi_g(i).ir_boundary = imgaussfilt(ir_b_vec,2);  
    master_struct_multi_g(i).dt_boundary = imgaussfilt(dt_b_vec,2);  
end

%% Plot range of achievable information rates
close all
alpha_index = master_struct_multi_bs(1).sweep_info_eq.a_index;
alpha_factor = master_struct_multi_bs(1).sweep_results_eq.rate_array(1,alpha_index);
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
        
% Same thing for decision times
topline = repelem(1e10,length(cw_index)-1);
close all
% make figure        
dt_space_fig_bs = figure;
hold on
feq = [];
for i = length(master_struct_multi_bs):-1:1
    cw_vec = master_struct_multi_bs(i).cw_boundary_eq;
    dt_vec = master_struct_multi_bs(i).dt_boundary_eq;
    
    fa = 0.75;
    ea = 0.75;
    if i == 1
        fa = 1;
    end
    if i == 1||i==5
        ea = 0.75;
    end
    cm = cmap_cmb(i,:);
   
    feq(end+1) = fill([cw_vec fliplr(cw_vec)],[dt_vec fliplr(topline)],...
                          cm,'FaceAlpha',fa,'EdgeAlpha',...
                          ea,'EdgeColor',brighten(cm,-0.5),'LineWidth',3); 
end

% plot IR noneq bound
% plot(master_struct_multi_bs(end).cw_boundary,master_struct_multi_bs(end).dt_boundary,':','Color',brighten(cmap_cmb(end,:),-0.25),'LineWidth',3)

% plot reference lines
yl = get(gca,'ylim');
yl(end) = 5e4;
yl(1) = 1;
plot(repmat(dt_cw_vec(plot_vec),2,1),repmat(yl',1,3),'-','Color','k','LineWidth',1)
plot(logspace(0,5),repelem(75,50),'-.','Color','k','LineWidth',1)
% errorbar(dt_cw_vec(plot_vec),dt_mean_vec(plot_vec),dt_mean_vec(plot_vec)-dt_lb_vec(plot_vec),dt_ub_vec(plot_vec)-dt_mean_vec(plot_vec),'o','LineWidth',1.5','Color','k')

legend(fliplr(feq),'n_b=1','n_b=2','n_b=3','n_b=4','n_b=5','Location','southeast')
xlabel('relative wrong factor concentration (w/c)');
ylabel('decision time (burst cycles)')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
set(gca,'yscale','log')
set(gca,'xscale','log')
set(gca,'xtick',[1 10 10^2 10^3 10^4 10^5])
xlim([1 1e5])
ylim([1 5e4])

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

dt_space_fig_bs.InvertHardcopy = 'off';
set(gcf,'color','w');
% ylim([0.25 1.75])

saveas(dt_space_fig_bs,[FigPath 'dt_vec_cw_bs.png'])
saveas(dt_space_fig_bs,[FigPath 'dt_vec_cw_bs.pdf']) 

%% Plot range of achievable information rates
topline = repelem(1e10,length(cw_index)-1);
close all

% make figure        
dt_space_fig_g = figure;
hold on
fneq = [];
for i = length(master_struct_multi_g)-2:-1:1
    cw_vec = master_struct_multi_g(i).cw_boundary;
    dt_vec = master_struct_multi_g(i).dt_boundary;
    
    fa = 0.75;
    ea = 0.75;
    if i == 1
        fa = 1;
        ea = 1;        
    end
    if i==5
        ea = 0.75;
    end
    cm = cmap_pu(2*i+1,:);
    if false
        cm = cmap2(2,:);
    end
    fneq(end+1) = fill([cw_vec fliplr(cw_vec)],[dt_vec fliplr(topline)],...
                          cm,'FaceAlpha',fa,'EdgeAlpha',...
                          ea,'EdgeColor',brighten(cm,-0.5),'LineWidth',3); 
end

plot(repmat(dt_cw_vec(plot_vec),2,1),repmat(yl',1,3),'-','Color','k','LineWidth',1)
plot(logspace(0,5),repelem(75,50),'-.','Color','k','LineWidth',1)

legend([fliplr(fneq)],'n_a=1','n_a=2','n_a=3','Location','northeast')
xlabel('relative wrong factor concentration (w/c)');
ylabel('decision time (burst cycles)')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
set(gca,'yscale','log')
set(gca,'xscale','log')
set(gca,'xtick',[1 10 10^2 10^3 10^4 10^5])
xlim([1 1e5])
ylim([1 5e4])

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

dt_space_fig_g.InvertHardcopy = 'off';
set(gcf,'color','w');
% ylim([0.25 1.75])

saveas(dt_space_fig_g,[FigPath 'dt_vec_cw_a.png'])
saveas(dt_space_fig_g,[FigPath 'dt_vec_cw_a.pdf'])

%% plot occupancy vs W/C
wc_axis = logspace(0,5);
% calculate occupancy vectors
spec_factor = 100;
spec_vec = [spec_factor ; spec_factor.^2 ; spec_factor.^3 ; spec_factor.^4];
occ_array = spec_vec ./ wc_axis;
plot_array = occ_array./(1+occ_array);

close all
occ_fig = figure;
hold on
for i = 1%:length(spec_vec)
    cm = cmap_pu(2*i,:);    
    plot(wc_axis,plot_array(i,:),'Color',brighten(cm,-0.1),'LineWidth',3.5)
end
    

xlabel('relative wrong factor concentration (w/c)');
ylabel('p_c/(p_c+p_w)')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
% set(gca,'yscale','log')
set(gca,'xscale','log')
set(gca,'xtick',[1 10 10^2 10^3 10^4 10^5])
xlim([1 1e5])
ylim([0 1.1])

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
occ_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(occ_fig,[FigPath 'occ_plot_eq.png'])
saveas(occ_fig,[FigPath 'occ_plot_eq.pdf'])

for i = 2:length(spec_vec)
    cm = cmap_pu(2*i,:);    
    plot(wc_axis,plot_array(i,:),'Color',brighten(cm,-0.6),'LineWidth',3.5)
end

legend('equilibrium','n_a=1','n_a=2','n_a=3','Location','southwest','Color','w')

saveas(occ_fig,[FigPath 'occ_plot.png'])
saveas(occ_fig,[FigPath 'occ_plot.pdf'])
%% Generate movie of transcription profile noise
mv_path = [FigPath 'mv_frames' filesep];
mkdir(mv_path);

c_vec = logspace(-2,2,1e3);
hill = 1;
tr_fun = c_vec.^hill./(c_vec.^hill + 1);

n_steps = 100;
sigma0 = 0.5;
noise_vec = normrnd(0,0.5,size(tr_fun));

close all
pu_dark = [190 189 216]/255;
pu_light = [229 221 238]/255;

for i = 1:n_steps
    temp_fig = figure('Visible','off');
    hold on
    sig_curr = 0.5/sqrt(i);
    n2 = normrnd(0,sig_curr/2,size(tr_fun));
    tr_noise = tr_fun+noise_vec/sqrt(i+1)+n2;
    tr_noise(tr_noise<0) = 0;
    plot(c_vec,tr_noise,'LineWidth',2,'Color',pu_light)
    plot(c_vec,tr_fun,'LineWidth',4,'Color',pu_dark)
    set(gca,'xscale','log')
    ylim([0 1.1])
    set(gca,'xtick',[],'ytick',[])
    saveas(temp_fig,[mv_path 'mv_frame_' sprintf('%03d',i) '.tif'])
end    

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
mdl_ir = fitlm(bs_vec'.^2, ir_vec_bs);
mdl_dt = fitlm(bs_vec(2:end)'.^-2, dt_vec_bs,'Intercept',false);

% extrapolate
bs_vec_pd = 0:100;
ir_vec_eq_pd = predict(mdl_ir,bs_vec_pd'.^2);
dt_vec_eq_pd = predict(mdl_dt,bs_vec_pd'.^-2);

% calculate dt coeficient
dtc = log((1-0.32)/.32)*(1-2*0.32);
% dt_vec_eq_pd = dtc./ir_vec_eq_pd;

figure;
hold on
scatter(bs_vec,ir_vec_bs);
plot(bs_vec_pd,ir_vec_eq_pd)
set(gca,'yscale','log')

figure;
hold on
scatter(bs_vec(2:end),dt_vec_bs);
plot(bs_vec_pd,dt_vec_eq_pd')
set(gca,'yscale','log')

