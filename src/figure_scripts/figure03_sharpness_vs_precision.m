% Plot results for sharpness vs precision for higher order models
clear 
close all
addpath(genpath('../utilities/'))

% %%%%%%%%%%%%%%%%  set relative read and write paths %%%%%%%%%%%%%%%%%%%%
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\';
DataPath = [DropboxFolder  'SweepOutput\sweeps02_sharpness_vs_precision' filesep ];
DataPathSupp = [DropboxFolder  'SweepOutput\sweeps01_info_vs_energy' filesep ];
FigPath = [DropboxFolder '\manuscript\sharpness_vs_precision' filesep];
mkdir(FigPath);

% %%%%%%%%%%%%%%%%  Set plot parameters and constants %%%%%%%%%%%%%%%%%%%%
n_plot = 3e3; % number of points to plot
n_samp = 5e4; % number of points to plot
markerAlpha = 0.5; % marker transparency
markerSize = 100; % marker size

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%% sharpness vs. precision for 4 state system %%%%%%%%%%%%%%%%
% set sweep options
sweep_options = {'n_sim',5,'n_seeds',5,'n_iters_max',50};

nStates = 4;
functionPath = ['../utilities/metricFunctions/n' num2str(nStates) '_OR/'];

% get metric names 
[~,~,metric_names_sym] = calculateMetricsSym_v2([]);
ir_index_sym = find(strcmp(metric_names_sym,'DecisionRateNorm'));
sharpness_index_sym = find(strcmp(metric_names_sym,'Sharpness'));
precision_index_sym = find(strcmp(metric_names_sym,'Precision'));
rate_index_sym = find(strcmp(metric_names_sym,'Production Rate'));
phi_index_sym = find(strcmp(metric_names_sym,'Phi'));
paramBounds = repmat([-5 ; 5],1,9);
paramBounds(:,1) = 0;
% run symbolic sweep
tic
[info_sym4, sweep_sym4] = ...
                param_sweep_multi_v3([sharpness_index_sym precision_index_sym],...
                                    functionPath,sweep_options{:},...
                                    'half_max_flag',false,'equilibrium_flag',false,...
                                    'TauCycleTime',1,'downsample_output',1,'paramBounds',paramBounds); 
                                  
[info_sym4_eq, sweep_sym4_eq] = ...
                param_sweep_multi_v3([sharpness_index_sym precision_index_sym],...
                                    functionPath,sweep_options{:},...
                                    'half_max_flag',false,'equilibrium_flag',true,...
                                    'TauCycleTime',1,'downsample_output',1,'paramBounds',paramBounds);                                   
toc

% calculate maximum info rates
%%%%%%%%%%
% neq
metric_array_sym_neq = vertcat(sweep_sym4.metric_array);   
sharpness_vec_neq = metric_array_sym_neq(:,sharpness_index_sym);
ir_vec_neq = metric_array_sym_neq(:,ir_index_sym);
precision_vec_neq = sqrt(exp(metric_array_sym_neq(:,precision_index_sym)));
plot_filter_neq = sharpness_vec_neq>=0;
plot_indices_neq = randsample(find(plot_filter_neq),n_plot,true);

%%%%%%%%%%%%%
% eq
metric_array_sym_eq = vertcat(sweep_sym4_eq.metric_array);   
sharpness_vec_eq = metric_array_sym_eq(:,sharpness_index_sym);
ir_vec_eq = metric_array_sym_eq(:,ir_index_sym);
precision_vec_eq = sqrt(exp(metric_array_sym_eq(:,precision_index_sym)));
plot_filter_eq = sharpness_vec_eq>=0;
plot_indices_eq = randsample(find(plot_filter_eq),n_plot,true);


% find best neq IR systems
rng(123);
n_ir = 100;
ir_thresh = 0.99*nanmax(ir_vec_neq);
% ir_99 = 0.99*nanmax(ir_vec_neq);
ir_options = find(ir_vec_neq>=ir_thresh&sharpness_vec_neq>=0);
ir_indices_plot = randsample(ir_options,n_ir,true);

% find most precise systems
p_99 = 0.99*nanmax(precision_vec_neq);
p_options = find(precision_vec_neq>=p_99&sharpness_vec_neq>=0);
p_indices_plot = randsample(p_options,n_ir,true);

% find sharpest systems
s_99 = 0.99*nanmax(sharpness_vec_neq);
s_options = find(sharpness_vec_neq>=s_99);
s_indices_plot = randsample(s_options,n_ir,true);

% ir_indices_plot = plot_indices_neq(ir_indices_plot);
% make figures
close 

sharp_v_prec_4_sc = figure;
cmap = brewermap(8,'Set2');
cmap_rd = brewermap(8,'Reds');

% plot eq and neq domains
hold on
scatter(sharpness_vec_neq(plot_indices_neq),precision_vec_neq(plot_indices_neq),markerSize,'MarkerFaceColor',...
            cmap(2,:),'MarkerEdgeColor','k','MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5)
scatter(sharpness_vec_eq(plot_indices_eq),precision_vec_eq(plot_indices_eq),markerSize,'MarkerFaceColor',...
            cmap(3,:),'MarkerEdgeColor','k','MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5) 
          
% plot position of 100 best IR gene circuits
scatter(sharpness_vec_neq(s_indices_plot),precision_vec_neq(s_indices_plot),...
        markerSize,'MarkerEdgeAlpha',1,'MarkerEdgeColor','k','MarkerFaceAlpha',0.5,...
      'MarkerFaceColor',brighten(cmap(5,:),-0.5)); 
    
scatter(sharpness_vec_neq(ir_indices_plot),precision_vec_neq(ir_indices_plot),...
        markerSize,'MarkerEdgeAlpha',1,'MarkerEdgeColor','k','MarkerFaceAlpha',0.5,...
      'MarkerFaceColor',brighten(cmap(8,:),-0.5)); 
    
scatter(sharpness_vec_neq(p_indices_plot),precision_vec_neq(p_indices_plot),...
    markerSize,'MarkerEdgeAlpha',1,'MarkerEdgeColor','k','MarkerFaceAlpha',0.5,...
  'MarkerFaceColor',cmap_rd(6,:)); 
       

xlabel('sharpness (S)');
ylabel('precision (P)')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

sharp_v_prec_4_sc.InvertHardcopy = 'off';
set(gcf,'color','w');

ylim([0 1.2])
xlim([0 2.2])

saveas(sharp_v_prec_4_sc,[FigPath 'sharp_vs_prec4_sc.png'])
% % saveas(sharp_v_prec_4_sc,[FigPath 'sharp_vs_prec4_sc.pdf'])    

%% Plot illustrative titration curves
% identify  systems near optimal sharpness for each gene circuit that
% attain maximum near HM
rate_bounds = [0.495 0.505];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find optimal neq systems
pd_vec_neq = metric_array_sym_neq(:,rate_index_sym);
rate_array_sym_neq = vertcat(sweep_sym4.rate_array);
pd_indices_neq = find(pd_vec_neq>=rate_bounds(1)&pd_vec_neq<=rate_bounds(2)&sharpness_vec_neq>=0);

% overall optimum
precision_ceiling = 0.4; % limit precision for illustrative purposes
[~,ir_i_neq] = nanmax(ir_vec_neq(pd_indices_neq).*(precision_vec_neq(pd_indices_neq)<=precision_ceiling));
ir_rates_neq = rate_array_sym_neq(pd_indices_neq(ir_i_neq),:);

% precision optimum
[p_max_neq,p_i_neq] = nanmax(precision_vec_neq(pd_indices_neq));
p_rates_neq = rate_array_sym_neq(pd_indices_neq(p_i_neq),:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find optimal eq system
pd_vec_eq = metric_array_sym_eq(:,rate_index_sym);
rate_array_sym_eq = vertcat(sweep_sym4_eq.rate_array);
pd_indices_eq = find(pd_vec_eq>=rate_bounds(1)&pd_vec_eq<=rate_bounds(2)&sharpness_vec_eq>=0);
[~,ir_i_eq] = nanmax(ir_vec_eq(pd_indices_eq));
ir_rates_eq = rate_array_sym_eq(pd_indices_eq(ir_i_eq),:);

% make arrays and calculate induction curves
c_vec = logspace(-2,2,1e3)';

ir_rate_array_neq = repmat(ir_rates_neq,length(c_vec),1);
ir_rate_array_neq(:,1) = c_vec;

p_rate_array_neq = repmat(p_rates_neq,length(c_vec),1);
p_rate_array_neq(:,1) = c_vec;

ir_rate_array_eq = repmat(ir_rates_eq,length(c_vec),1);
ir_rate_array_eq(:,1) = c_vec;

% map to correct subfolder
functionPath = info_sym4.functionPath;
addpath(genpath(functionPath));

% calculate production rates and noise profiles
paramCellNeq = mat2cell(ir_rate_array_neq,size(ir_rate_array_neq,1),ones(1,size(ir_rate_array_neq,2)));
paramCellNeqP = mat2cell(p_rate_array_neq,size(p_rate_array_neq,1),ones(1,size(p_rate_array_neq,2)));
paramCellEq = mat2cell(ir_rate_array_eq,size(ir_rate_array_eq,1),ones(1,size(ir_rate_array_eq,2)));

% neq
r_curve_neq = productionRateFunction(paramCellNeq{:});
p_curve_neq = sqrt(intrinsicVarianceFunction(paramCellNeq{:}));
% prec neq
r_curve_neq_p = productionRateFunction(paramCellNeqP{:});
p_curve_neq_p = sqrt(intrinsicVarianceFunction(paramCellNeqP{:}));
% eq
r_curve_eq = productionRateFunction(paramCellEq{:});
p_curve_eq = sqrt(intrinsicVarianceFunction(paramCellEq{:}));

% Make plot 
close all
c1 = info_sym4.cr1;
c0 = info_sym4.cr0;

induction_plots_sym = figure;
cmap = brewermap(8,'Set2');
cmap_rd = brewermap(8,'Reds');

hold on
n_cycles = 25;%60/5;
% calculate upper and lower bounds 
ubeq = r_curve_eq + p_curve_eq/sqrt(n_cycles);
lbeq = r_curve_eq - p_curve_eq/sqrt(n_cycles);

ubneq = r_curve_neq + p_curve_neq/sqrt(n_cycles);
lbneq = r_curve_neq - p_curve_neq/sqrt(n_cycles);

ubp = r_curve_neq_p + p_curve_neq_p/sqrt(n_cycles);
lbp = r_curve_neq_p - p_curve_neq_p/sqrt(n_cycles);

% plot verticle lines denoting target concentrations
plot(repelem(c0,100),linspace(0,1,100),'-.k','LineWidth',1.5)
plot(repelem(c1,100),linspace(0,1,100),'-.k','LineWidth',1.5)

% make area plots to indicate error range for each
ceq = cmap(3,:);
fill([c_vec' fliplr(c_vec')],[lbeq' fliplr(ubeq')],ceq,'FaceAlpha',.1,'EdgeAlpha',0.75,'EdgeColor',ceq);
peq = plot(c_vec,r_curve_eq,'Color',cmap(3,:),'LineWidth',2);

cneq = brighten(cmap(5,:),-0.5);
fill([c_vec' fliplr(c_vec')],[ubneq' fliplr(lbneq')],cneq,'FaceAlpha',.1,'EdgeAlpha',0.75,'EdgeColor',cneq);
pneq = plot(c_vec,r_curve_neq,'Color',brighten(cmap(5,:),-0.5),'LineWidth',2);

cp = cmap_rd(6,:);
fill([c_vec' fliplr(c_vec')],[ubp' fliplr(lbp')],cp,'FaceAlpha',.1,'EdgeAlpha',0.75,'EdgeColor',cp);
pp = plot(c_vec,r_curve_neq_p,'Color',cp,'LineWidth',2);



set(gca,'xscale','log')

% set(gca,'xtick',[0.2 1 5])
xlim([0.1 10])
ylim([0 1])
xlabel('activator concentration (c)');
ylabel('transcription rate (r)')

% legend([peq pneq pp],'equilibrium','non-eq. (sharp)','non-eq. (precise)','Location','northwest')
grid on

set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

induction_plots_sym.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(induction_plots_sym,[FigPath 'motif_induction_4state.png'])
saveas(induction_plots_sym,[FigPath 'motif_induction_4state.pdf'])    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%% Examine impact of adding binding sites %%%%%%%%%%%%%%%%%
% get metric names for numeric sweeps
[~,~,metric_names_num] = calculateMetricsNumeric_v3([]);

ir_index_num = find(strcmp(metric_names_num,'IR'));
sharpness_index_num = find(strcmp(metric_names_num,'Sharpness'));
precision_index_num = find(strcmp(metric_names_num,'Precision'));
rate_index_num = find(strcmp(metric_names_num,'ProductionRate'));
phi_index_num = find(strcmp(metric_names_num,'Phi'));
tau_index_num = find(strcmp(metric_names_num,'TauCycle'));


% get list of sweep results files with only 1 genera TF reaction
multi_bs_sweep_files = dir([DataPath 'sweep_results*g01*neq*']);
multi_bs_info_files = dir([DataPath 'sweep_info*g01*neq*']);

% load
master_struct_multi_bs = struct;
for f = 1:length(multi_bs_sweep_files)
  
    load([DataPath multi_bs_sweep_files(f).name])
    load([DataPath multi_bs_info_files(f).name])
    
    master_struct_multi_bs(f).sweep_results = sim_results;
    master_struct_multi_bs(f).sweep_info = sim_info;
end    

% calculate sharpness vs precision boundaries 
% use precision and sharpness to eliminate outliers due to numerical 
% precision for now. Need to revisit this
% p_max = [2.25 2.25 2.25 2.25 1.62.^2];
s_min = 0.005;

rng(321);

for i = 1:length(master_struct_multi_bs)
    % extract vectors
    metric_array = master_struct_multi_bs(i).sweep_results.metric_array;
    rate_array = master_struct_multi_bs(i).sweep_results.rate_array;
    sharpness_vec = metric_array(:,sharpness_index_num).^2;
    precision_vec = metric_array(:,precision_index_num).^2;
    ir_vec = metric_array(:,ir_index_num);

    if i > 4
        ft2 = sqrt(sharpness_vec)<=s_min;% NL: screen out a handful of outliers with very small sharpness and anamaolously high precision
    else
        ft2 = false(size(sharpness_vec));
    end
    options = find(~ft2);
    
    % identify boundary points
    b_points = boundary(sqrt(sharpness_vec(options)),sqrt(precision_vec(options)),0.9);
    
    % store
    master_struct_multi_bs(i).bound_points = options(b_points);
    master_struct_multi_bs(i).sharpness_boundary = sharpness_vec(options(b_points));
    master_struct_multi_bs(i).precision_boundary = precision_vec(options(b_points));
    
    % select N top performers to plot
    ir_thresh = prctile(ir_vec(options),99.9);%0.98*nanmax(ir_vec(options));
    ir_options = find(ir_vec(options)>=ir_thresh);
    ns = min([n_ir length(options(ir_options))]);
    if length(ir_options) > 1
        ir_indices = randsample(options(ir_options),ns,true);
    else
        ir_indices = options(ir_options);
    end
    master_struct_multi_bs(i).sharp_ir = sharpness_vec(ir_indices);
    master_struct_multi_bs(i).prec_ir = precision_vec(ir_indices);
    
    % save max sharpness and precision
    master_struct_multi_bs(i).sharp_max = nanmax(sharpness_vec(options));
    master_struct_multi_bs(i).prec_max = nanmax(precision_vec(options));
end


close all

% set plot parameters
% markerSize = 50;
rng(231);
alphaFactor = 0.25;

% Define colormaps for use throughout
cmap_pu = brewermap(8,'Purples');
cmap_rd = brewermap(8,'Reds');
cmap_bu = brewermap(8,'Blues');
cmap_gre = brewermap(8,'Greens');
cmap_gra = brewermap(8,'Greys');

close all
color_ind = 5;
cmap = [cmap_pu(3,:); cmap_gre(color_ind,:); cmap_rd(color_ind,:); ...
          cmap_bu(color_ind,:) ; cmap_gra(color_ind,:)];
       
% make figure        
sharp_vs_prec_bs = figure;
hold on
% plot area vectors
for i = length(master_struct_multi_bs):-1:1                      
    % downsample vectors    
    s_vec = sqrt(master_struct_multi_bs(i).sharpness_boundary);
    p_vec = sqrt(master_struct_multi_bs(i).precision_boundary);
    fill([s_vec fliplr(s_vec)],[p_vec fliplr(p_vec)],...
                                        cmap(i,:),'FaceAlpha',alphaFactor*2,'EdgeAlpha',...
                                        alphaFactor*2,'EdgeColor',brighten(cmap(i,:),-0.5),'LineWidth',2);    
end

for i = length(master_struct_multi_bs):-1:1  
    s_vec = sqrt(master_struct_multi_bs(i).sharp_ir);
    p_vec = sqrt(master_struct_multi_bs(i).prec_ir);
    scatter(s_vec,p_vec,markerSize,'MarkerEdgeAlpha',1,'MarkerEdgeColor','k','MarkerFaceAlpha',0.75,...
                  'MarkerFaceColor',cmap(i,:)); 
end    

xlabel('sharpness (S)');
ylabel('precision (P)')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

sharp_vs_prec_bs.InvertHardcopy = 'off';
set(gcf,'color','w');
ylim([0.25 1.75])

saveas(sharp_vs_prec_bs,[FigPath 'sharp_vs_prec_bs.png'])
saveas(sharp_vs_prec_bs,[FigPath 'sharp_vs_prec_bs.pdf']) 
%%

% make norm figure        
eq_s_vec = (1:5);
eq_p_vec = sqrt(1/2);
close all
sharp_vs_prec_bs = figure;
hold on
% plot area vectors
for i = length(master_struct_multi_bs):-1:1                     
    % downsample vectors    
    s_vec = sqrt(master_struct_multi_bs(i).sharpness_boundary)/eq_s_vec(i);
    p_vec = sqrt(master_struct_multi_bs(i).precision_boundary)/eq_p_vec;
    fa = 0.05;
    if i == 1
        fa = 0.2;
    end
    fill([s_vec fliplr(s_vec)],[p_vec fliplr(p_vec)],...
                                        cmap(i,:),'FaceAlpha',fa,'EdgeAlpha',...
                                        alphaFactor*2,'EdgeColor',brighten(cmap(i,:),-0.5),'LineWidth',3);    
end

for i = 1:length(master_struct_multi_bs) 
    s_vec = sqrt(master_struct_multi_bs(i).sharp_ir)/eq_s_vec(i);
    p_vec = sqrt(master_struct_multi_bs(i).prec_ir)/eq_p_vec;
    scatter(s_vec,p_vec,markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',0.25,...
                  'MarkerFaceColor',cmap(i,:)); 
end    


xlabel('non-equilibrium sharpness gain (S/S_{eq})');
ylabel('non-equilibrium precision gain  (P/P_{eq})')
xlim([0 2.1])
ylim([0 2.5])
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

sharp_vs_prec_bs.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(sharp_vs_prec_bs,[FigPath 'sharp_vs_prec_norm_bs.png'])
saveas(sharp_vs_prec_bs,[FigPath 'sharp_vs_prec_norm_bs.pdf']) 

% Plot sharpness and precision as a function of nbs

xlbs = [0.5 5.5];

close all
sharp_bs = figure;

hold on

% make plots
n_plot = 500;
plot(linspace(xlbs(1), xlbs(2)),2*linspace(xlbs(1), xlbs(2)),':','Color','k','LineWidth',2)
for i = 1:length(master_struct_multi_bs)
    s_vec = master_struct_multi_bs(i).sweep_results.metric_array(:,sharpness_index_num);
    s_options = find(s_vec>=0);
    plot_ids = randsample(s_options,n_plot,true,sqrt(s_vec));
    g_temp = i + normrnd(0,0.025,n_plot,1);
    scatter(g_temp,s_vec(s_options(plot_ids)),25,'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.1);
end

ylabel('normalized sharpness (S)')
% ylim([2 5])
xlim(xlbs)
xlabel('number of binding sites');

grid on
box on
set(gca,'xtick',1:4)
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
% ylim([0 0.4])
sharp_bs.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(sharp_bs,[FigPath 'sharp_vs_bs_sc.png'])
saveas(sharp_bs,[FigPath 'sharp_vs_bs_sc.pdf']) 

% note that I'm using low Phi as a proxy for eq performance
s_max_vec = sqrt([master_struct_multi_bs.sharp_max]);
p_max_vec = sqrt([master_struct_multi_bs.prec_max]);
bs_pd_vec = 1:5;
xl = [0.5 5.5];
close all
sharp_bs = figure;

hold on

% make plots
plot(linspace(xl(1),xl(2)),2*linspace(xl(1),xl(2)),'--','Color','k','LineWidth',2)
scatter(bs_pd_vec,s_max_vec,75,'MarkerFaceColor',cmap_gre(4,:),'MarkerEdgeColor','k');

ylabel('maximum achivable sharpness (S)')
% ylim([2 10])

xlabel('number of binding sites');

grid on
box on

set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.XAxis(1).Color = 'k';
ax.YAxis(1).Color = 'k';
xlim(xl);
sharp_bs.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(sharp_bs,[FigPath 'sharp_vs_bs.png'])
saveas(sharp_bs,[FigPath 'sharp_vs_bs.pdf']) 

%%
prec_bs = figure;

hold on

plot(bs_pd_vec,p_max_vec,'-','Color','k','LineWidth',2)
scatter(bs_pd_vec,p_max_vec,75,'s','MarkerFaceColor',cmap_rd(4,:),'MarkerEdgeColor','k');

ylabel('maximum achivable precision (P)')
ylim([0.5 2])

xlabel('number of binding sites');

grid on
box on

set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
% ylim([0 0.4])
sharp_bs.InvertHardcopy = 'off';
set(gcf,'color','w');
xlim(xl)
saveas(sharp_bs,[FigPath 'prec_vs_bs.png'])
saveas(sharp_bs,[FigPath 'prec_vs_bs.pdf']) 


%% Plot illustrative induction curves for 
DataPath2 = [DropboxFolder  'SweepOutput\sweeps02B_sharpness_vs_rate' filesep];

multi_bs_rate_files = dir([DataPath2 'sweep_results*g01*neq*']);

% identify  systems near optimal sharpness for each gene circuit that
% attain maximum near HM
rng(123);
for i = 1:length(multi_bs_rate_files)
  
    % read file into working space
    load([DataPath2 multi_bs_rate_files(i).name])
    master_struct_multi_bs(i).sweep_results_rate = sim_results;
       
    s_vec = master_struct_multi_bs(i).sweep_results_rate.metric_array(:,sharpness_index_num);
    r_vec = master_struct_multi_bs(i).sweep_results_rate.metric_array(:,rate_index_num);
    phi_vec = master_struct_multi_bs(i).sweep_results_rate.metric_array(:,phi_index_num);
    r_bounds = [0.49 0.51];
    
    options = find(r_vec>=r_bounds(1) & r_vec<=r_bounds(2));
    [~,sharp_i] = nanmax(s_vec(options));
    master_struct_multi_bs(i).sharp_rates = master_struct_multi_bs(i).sweep_results_rate.rate_array(options(sharp_i),:);        
end    

% generate predicted response curve
c_vec = logspace(-1, 1,100);
numerical_precision = 5;
induction_array_bs = NaN(length(c_vec),length(master_struct_multi_bs));

for i = 1:length(master_struct_multi_bs)
  
    % add correct function to working path
    sim_info = master_struct_multi_bs(i).sweep_info;
    functionPath = sim_info.functionPath;
    slashes = strfind(functionPath,'\');
    simName = functionPath(slashes(end-1)+1:slashes(end)-1);
    rmpath(genpath('../utilities/metricFunctions/'));
    addpath(genpath(['../utilities/metricFunctions/numeric/' simName]));
    
    % extract rates
    param_vec = master_struct_multi_bs(i).sharp_rates;
    
    for c = 1:length(c_vec)
        param_temp = param_vec;
        param_temp(1) = c_vec(c);
        valCellCS = mat2cell(param_temp,size(param_temp,1),ones(1,size(param_temp,2)));   
        
        Q_num = RSymFun(valCellCS{:});

        ss_short = calculate_ss_num(Q_num,numerical_precision);  
        induction_array_bs(c,i) = sum(ss_short(sim_info.activeStateFilter));
    end
end

c1 = master_struct_multi_bs(1).sweep_info.cr1;
c0 = master_struct_multi_bs(1).sweep_info.cr0;

% Make plot 
sharpness_plots_bs = figure;
hold on
% plot area vectors
p = [];
for i = 1:length(master_struct_multi_bs)                     
    % downsample vectors    
    p(i) = plot(c_vec,induction_array_bs(:,i),'Color',cmap(i,:),'LineWidth',3);
   
end

% plot verticle lines denoting target concentrations
plot(repelem(c0,100),linspace(0,1,100),':k','LineWidth',2)
plot(repelem(c1,100),linspace(0,1,100),':k','LineWidth',2)

legend(p,'1 bs','2 bs','3 bs','4 bs','5 bs','Location','southeast');
xlim([0.1 10])

set(gca,'xscale','log')

xlabel('activator concentration (c)');
ylabel('transcription rate (r)')

grid on
% set(gca,'xtick',[0.2 1 5])
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

sharpness_plots_bs.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(sharpness_plots_bs,[FigPath 'sharp_curves_bs.png'])
saveas(sharpness_plots_bs,[FigPath 'sharp_curves_bs.pdf']) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%% Look at impact of multiple molecular steps %%%%%%%%%%%%%%%%%%%%

% get list of sweep results files with only 1 genera TF reaction
multi_g_sweep_files = dir([DataPath 'sweep_results_s01_ns00_g0*neq*']);
multi_g_info_files = dir([DataPath 'sweep_info_s01_ns00_g0*neq*']);
  
% load
master_struct_multi_g = struct;
for f = 1:length(multi_g_sweep_files)
  
    load([DataPath multi_g_sweep_files(f).name])
    load([DataPath multi_g_info_files(f).name])
    
    master_struct_multi_g(f).sweep_results = sim_results;
    master_struct_multi_g(f).sweep_info = sim_info;
end    

% Supplement with Phi vs. IR sweep info
% This is used to screen outliers from S vs. P sweep and to fill out
% under-represented regions
multi_g_sweep_files_supp = dir([DataPathSupp 'sweep_results_s01_ns00_g0*']);
multi_g_info_files_supp = dir([DataPathSupp 'sweep_info_s01_ns00_g0*']);
for f = 1:length(multi_g_sweep_files)
  
    load([DataPathSupp multi_g_sweep_files_supp(f).name])
    load([DataPathSupp multi_g_sweep_files_supp(f).name])
    
    master_struct_multi_g(f).sweep_results_supp = sim_results;
    master_struct_multi_g(f).sweep_info_supp = sim_info;
end    
% calculate sharpness vs precision boundaries 
% use IR to reject few outliers for now
% n_ir = 50;
% p_max = 1.6;
% s_max = 6^2;
% ir_max_vec = [0.015 0.03 0.045 0.06 0.075];
rng(321);

for i = 1:length(master_struct_multi_g)-1
    % extract vectors
    metric_array = master_struct_multi_g(i).sweep_results.metric_array; 
    metric_array_supp = master_struct_multi_g(i).sweep_results_supp.metric_array; 
    sharpness_vec = vertcat(metric_array(:,sharpness_index_num).^2,metric_array_supp(:,sharpness_index_num).^2);    
    precision_vec = vertcat(metric_array(:,precision_index_num).^2,metric_array_supp(:,precision_index_num).^2);
    
    ir_vec = vertcat(metric_array(:,ir_index_num),metric_array_supp(:,ir_index_num))*log2(exp(1));
    ir_ceil = nanmax(metric_array_supp(:,ir_index_num))*log2(exp(1)); % anything higher than this is a num precision error
    options = 1:length(sharpness_vec);%find(sharpness_vec<=s_max&precision_vec<=p_max&ir_vec<=ir_ceil);
    
    % identify boundary points
    if i > 1
        b_points = boundary(sqrt(sharpness_vec(options)),sqrt(precision_vec(options)),0.6);
    else
        b_points = boundary(sqrt(sharpness_vec(options)),sqrt(precision_vec(options)),0.3);
    end
    % store
    master_struct_multi_g(i).bound_points = options(b_points);
    master_struct_multi_g(i).sharpness_boundary = sharpness_vec(options(b_points));
    master_struct_multi_g(i).precision_boundary = precision_vec(options(b_points));
    
    % select N top performers to plot
    ir_thresh = prctile(ir_vec(options),99.9);%0.99*nanmax(ir_vec(options));
    ir_options = find(ir_vec(options)>=ir_thresh);

    if length(ir_options)>1
        ir_indices = randsample(options(ir_options),min([n_ir,length(ir_options)]),false);
    else
        ir_indices = options(ir_options);
    end
    master_struct_multi_g(i).sharp_ir = sharpness_vec(ir_indices);
    master_struct_multi_g(i).prec_ir = precision_vec(ir_indices);
    
    % save max sharpness and precision
    master_struct_multi_g(i).sharp_max = nanmax(sharpness_vec(options));
    master_struct_multi_g(i).prec_max = nanmax(precision_vec(options));
    
end

% make figure        
close all
sharp_vs_prec_g = figure;
hold on
% plot area vectors
for i = length(master_struct_multi_g)-1:-1:1                      
    % downsample vectors    
    s_vec = sqrt(master_struct_multi_g(i).sharpness_boundary);
    p_vec = sqrt(master_struct_multi_g(i).precision_boundary);
    fill([s_vec fliplr(s_vec)],[p_vec fliplr(p_vec)],...
                                        cmap_pu(2+i,:),'FaceAlpha',.05,'EdgeAlpha',...
                                        1,'EdgeColor',brighten(cmap_pu(2+i,:),-0.25),'LineWidth',3);    
end

for i = length(master_struct_multi_g)-1:-1:1  
    s_vec = sqrt(master_struct_multi_g(i).sharp_ir);
    p_vec = sqrt(master_struct_multi_g(i).prec_ir);
    scatter(s_vec,p_vec,markerSize,'MarkerEdgeAlpha',1,'MarkerEdgeColor','k','MarkerFaceAlpha',1,...
                  'MarkerFaceColor',cmap_pu(2+i,:)); 
end    

xlabel('sharpness (S)');
ylabel('precision (P)')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

sharp_vs_prec_g.InvertHardcopy = 'off';
set(gcf,'color','w');
ylim([0 1.6])

saveas(sharp_vs_prec_g,[FigPath 'sharp_vs_prec_g.png'])
saveas(sharp_vs_prec_g,[FigPath 'sharp_vs_prec_g.pdf']) 


% make norm figure        
eq_s_vec = 1;
markerSize = 100;
eq_p_vec = sqrt(1/2);

sharp_vs_prec_g = figure;
hold on
% plot area vectors
for i = length(master_struct_multi_g)-1:-1:1                    
    % downsample vectors    
    s_vec = sqrt(master_struct_multi_g(i).sharpness_boundary)/eq_s_vec;
    p_vec = sqrt(master_struct_multi_g(i).precision_boundary)/eq_p_vec;
    if i > 1
      fa = 0.05;
    else 
      fa = 0.2;
    end
    fill([s_vec fliplr(s_vec)],[p_vec fliplr(p_vec)],...
                                        cmap_pu(i+2,:),'FaceAlpha',fa,'EdgeAlpha',...
                                        1,'EdgeColor',brighten(cmap_pu(i+2,:),-0.5),'LineWidth',3);    
end

for i = length(master_struct_multi_g)-1:-1:1  
    s_vec = sqrt(master_struct_multi_g(i).sharp_ir)/eq_s_vec;
    p_vec = sqrt(master_struct_multi_g(i).prec_ir)/eq_p_vec;
    scatter(s_vec,p_vec,markerSize,'MarkerEdgeAlpha',1,'MarkerEdgeColor','k','MarkerFaceAlpha',0.75,...
                  'MarkerFaceColor',cmap_pu(2+i,:)); 
end    

xlabel('non-equilibrium sharpness gain (S/S_{eq})');
ylabel('non-equilibrium precision gain (P/P_{eq})')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
ylim([0 2])
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

sharp_vs_prec_g.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(sharp_vs_prec_g,[FigPath 'sharp_vs_prec_norm_g.png'])
saveas(sharp_vs_prec_g,[FigPath 'sharp_vs_prec_norm_g.pdf']) 

% Plot max sharpness and precision
% note that I'm using low Phi as a proxy for eq performance

xlg = [0.5 4.5];

close all
sharp_g = figure;

hold on

% make plots
n_plot = 500;
plot(linspace(xlg(1), xlg(2)),1+linspace(xlg(1), xlg(2)),':','Color','k','LineWidth',2)
for i = 1:length(master_struct_multi_g)-1
    s_vec = master_struct_multi_g(i).sweep_results.metric_array(:,sharpness_index_num);
    s_options = find(s_vec>=0);
    plot_ids = randsample(s_options,n_plot,true,sqrt(s_vec));
    g_temp = i + normrnd(0,0.025,n_plot,1);
    scatter(g_temp,s_vec(s_options(plot_ids)),25,'MarkerFaceColor',cmap_pu(i*2,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.1);
end

ylabel('sharpness (S)')
% ylim([2 5])
xlim(xlg)
xlabel('number of locus conformations');

grid on
box on
set(gca,'xtick',1:4)
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
% ylim([0 0.4])
sharp_g.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(sharp_g,[FigPath 'sharp_vs_g_sc.png'])
saveas(sharp_g,[FigPath 'sharp_vs_g_sc.pdf']) 

% note that I'm using low Phi as a proxy for eq performance
s_max_vec_g = sqrt([master_struct_multi_g(1:4).sharp_max]);
p_max_vec_g = sqrt([master_struct_multi_g(1:4).prec_max]);
g_vec = 1:4;
xlg = [0.5 4.5];

close all
sharp_g = figure;

hold on

% make plots
plot(linspace(xlg(1), xlg(2)),1+linspace(xlg(1), xlg(2)),'.','Color','k','LineWidth',2)
scatter(g_vec,s_max_vec_g,75,'MarkerFaceColor',cmap_gre(4,:),'MarkerEdgeColor','k');

ylabel('maximum  sharpness (S)')
% ylim([2 5])
xlim(xlg)
xlabel('number of locus conformations');

grid on
box on
set(gca,'xtick',1:4)
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
% ylim([0 0.4])
sharp_g.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(sharp_g,[FigPath 'sharp_vs_g.png'])
saveas(sharp_g,[FigPath 'sharp_vs_g.pdf']) 

%
prec_g = figure;

hold on

% make plots

plot(g_vec,p_max_vec_g,'-','Color','k','LineWidth',2)
scatter(g_vec,p_max_vec_g,75,'s','MarkerFaceColor',cmap_rd(4,:),'MarkerEdgeColor','k');

ylabel('maximum precision (P)')
ylim([.75 1.4])
xlim(xlg)

xlabel('number of general TF reactions');

grid on
box on
set(gca,'xtick',1:4)
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
% ylim([0 0.4])
prec_g.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(prec_g,[FigPath 'prec_vs_g.png'])
saveas(prec_g,[FigPath 'prec_vs_g.pdf']) 


%% Plot illustrative induction curves for different Ng

multi_g_rate_files = dir([DataPath2 'sweep_results_s01_ns00_g0*neq*']);


% identify  systems near optimal sharpness for each gene circuit that
% attain maximum near HM
rng(123);
for i = 1:length(master_struct_multi_g)

    r_bounds = [0.49 0.51];
    
    % read file into working space
    load([DataPath2 multi_g_rate_files(i).name])
    master_struct_multi_g(i).sweep_results_rate = sim_results;
       
    s_vec = master_struct_multi_g(i).sweep_results_rate.metric_array(:,sharpness_index_num);
    r_vec = master_struct_multi_g(i).sweep_results_rate.metric_array(:,rate_index_num);
    
    s_99 = 0.95*nanmax(s_vec);
    options = find(s_vec>=s_99 & r_vec>=r_bounds(1) & r_vec<=r_bounds(2));

    if isempty(options)
        error('no options')        
    elseif length(options) == 1
        master_struct_multi_g(i).sharp_rates = master_struct_multi_g(i).sweep_results_rate.rate_array(options,:);
    else
        ind = randsample(options,1);
        master_struct_multi_g(i).sharp_rates = master_struct_multi_g(i).sweep_results_rate.rate_array(ind,:);
    end
end    

% generate predicted response curve
c_vec = logspace(-1, 1,100);
numerical_precision = 5;
induction_array_g = NaN(length(c_vec),length(master_struct_multi_g));

for i = 1:length(master_struct_multi_g)
  
    % add correct function to working path
    sim_info = master_struct_multi_g(i).sweep_info;
    functionPath = sim_info.functionPath;
    slashes = strfind(functionPath,'\');
    simName = functionPath(slashes(end-1)+1:slashes(end)-1);
    rmpath(genpath('../utilities/metricFunctions/'));
    addpath(genpath(['../utilities/metricFunctions/numeric/' simName]));
    
    % extract rates
    param_vec = master_struct_multi_g(i).sharp_rates;
    
    for c = 1:length(c_vec)
        param_temp = param_vec;
        param_temp(1) = c_vec(c);
        valCellCS = mat2cell(param_temp,size(param_temp,1),ones(1,size(param_temp,2)));   
        
        Q_num = RSymFun(valCellCS{:});

        ss_short = calculate_ss_num(Q_num,numerical_precision);  
        induction_array_g(c,i) = sum(ss_short(sim_info.activeStateFilter));
    end
end

c1 = master_struct_multi_g(1).sweep_info.cr1;
c0 = master_struct_multi_g(1).sweep_info.cr0;

% Make plot 
sharpness_plots_g = figure;
hold on
% plot area vectors
p = [];
for i = 1:length(master_struct_multi_g)-1                     
    % downsample vectors    
    p(i) = plot(c_vec,induction_array_g(:,i),'Color',brighten(cmap_pu(2+i,:),-0.25),'LineWidth',3);
   
end

% plot verticle lines denoting target concentrations
plot(repelem(c0,100),linspace(0,1,100),':k','LineWidth',1.5)
plot(repelem(c1,100),linspace(0,1,100),':k','LineWidth',1.5)

legend(p,'1 LC','2 LC','3 LC','4 LC','Location','southeast');
xlim([0.1 10])

set(gca,'xscale','log')
xlabel('activator concentration (c)');
ylabel('transcription rate (r)')

% grid on

set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

sharpness_plots_g.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(sharpness_plots_g,[FigPath 'sharp_curves_g.png'])
saveas(sharpness_plots_g,[FigPath 'sharp_curves_g.pdf']) 

%%
DataPath3 = [DropboxFolder  'SweepOutput\sweeps02_sharpness_vs_precision' filesep];

multi_bs_rate_files = dir([DataPath3 'sweep_results*g01*_eq*']);
multi_bs_info_files = dir([DataPath3 'sweep_info*g01*_eq*']);
nb_ind = 3;
% identify  systems near optimal sharpness for each gene circuit that
% attain maximum near HM
rng(123);

% read file into working space
clear sim_results
clear sim_info
load([DataPath3 multi_bs_rate_files(nb_ind).name])
load([DataPath3 multi_bs_info_files(nb_ind).name])
master_struct_multi_bs(nb_ind).sweep_results_rate = sim_results;

s_vec = sim_results.metric_array(:,sharpness_index_num);
r_vec = sim_results.metric_array(:,rate_index_num);
r_bounds = [0.49 0.51];

options = find(r_vec>=r_bounds(1) & r_vec<=r_bounds(2));
[~,sharp_i] = nanmax(s_vec(options));
sharp_rates_bs = sim_results.rate_array(options(sharp_i),:);        

% get predicted induction curve
% add correct function to working path
functionPath = sim_info.functionPath;
slashes = strfind(functionPath,'\');
simName = functionPath(slashes(end-1)+1:slashes(end)-1);
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(['../utilities/metricFunctions/numeric/' simName]));

% extract rates
param_vec = sharp_rates_bs;
bs_pd_vec = NaN(size(c_vec));

for c = 1:length(c_vec)
    param_temp = param_vec;
    param_temp(1) = c_vec(c);
    valCellCS = mat2cell(param_temp,size(param_temp,1),ones(1,size(param_temp,2)));   

    Q_num = RSymFun(valCellCS{:});

    ss_short = calculate_ss_num(Q_num,numerical_precision);  
    bs_pd_vec(c) = sum(ss_short(sim_info.activeStateFilter));
end

%% Make figure comparing multi-bs and multi-step
close all

na_ind = 2;


% Make plot 
sharpness_plots_na_vs_nb = figure;
hold on
% plot area vectors
p = [];
              
% downsample vectors    
p_na = plot(c_vec,induction_array_g(:,na_ind),'Color',cmap_pu(2+na_ind,:),'LineWidth',3);
p_nb = scatter(c_vec(1:3:end),bs_pd_vec(1:3:end),'o','MarkerFaceColor',cmap(nb_ind,:),'MarkerEdgeColor','k');


% plot verticle lines denoting target concentrations
% plot(repelem(c0,100),linspace(0,1,100),':k','LineWidth',1.5)
% plot(repelem(c1,100),linspace(0,1,100),':k','LineWidth',1.5)

legend([p_na p_nb],'N_a= 3 (non-eq.)','N_b= 3 (eq.)','Location','northwest');
xlim([0.1 10])

set(gca,'xscale','log')
xlabel('activator concentration (c)');
ylabel('transcription rate (r)')

% grid on
box on
set(gca,'FontSize',24)
set(gca,'Color',[228,221,209]/255) 
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
ax.LineWidth = 3;
sharpness_plots_na_vs_nb.InvertHardcopy = 'off';
set(gcf,'color','w');
% 
saveas(sharpness_plots_na_vs_nb,[FigPath 'na_vs_nb.png'])
saveas(sharpness_plots_na_vs_nb,[FigPath 'na_vs_nb.pdf']) 