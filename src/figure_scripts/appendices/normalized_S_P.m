% Make figures illustrating normalized sharpness and precision metrics
clear 
close all
addpath(genpath('../utilities/'))

% %%%%%%%%%%%%%%%%  set relative read and write paths %%%%%%%%%%%%%%%%%%%%
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\';
FigPath = [DropboxFolder '\manuscript\Appendices' filesep];
mkdir(FigPath);

% %%%%%%%%%%%%%%%%  Set plot parameters and constants %%%%%%%%%%%%%%%%%%%%
% n_plot = 3e3; % number of points to plot
n_samp = 5e4; % number of points to plot
markerAlpha = 0.5; % marker transparency
markerSize = 75; % marker size

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%% sharpness vs. precision for 4 state system %%%%%%%%%%%%%%%%
% set sweep options
sweep_options = {'n_sim',1,'n_seeds',5,'n_iters_max',50};

nStates = 4;
functionPath = ['../../utilities/metricFunctions/n' num2str(nStates) '_OR/'];

% get metric names 
[~,~,metric_names_sym] = calculateMetricsSym_v2([]);
ir_index_sym = find(strcmp(metric_names_sym,'DecisionRateNorm'));
sharpness_index_sym = find(strcmp(metric_names_sym,'Sharpness'));
precision_index_sym = find(strcmp(metric_names_sym,'Precision'));
sharpness_raw_index_sym = find(strcmp(metric_names_sym,'SharpnessRaw'));
precision_raw_index_sym = find(strcmp(metric_names_sym,'PrecisionRaw'));
rate_index_sym = find(strcmp(metric_names_sym,'Production Rate'));
paramBounds = repmat([-5 ; 5],1,9);
paramBounds(:,1) = 0;

% initialize structure to store results
sweep_struct_master = struct;
sweep_id_cell = {'p','NP','s','NS'};
y_axis_cell = {'precision (p)','normalized precision (P)','sharpness (s)','normalized sharpness (S)'};
sweep_index_vec = [precision_raw_index_sym precision_index_sym sharpness_raw_index_sym sharpness_index_sym];

for s = 1:length(sweep_index_vec)
    % run symbolic sweep
    tic
    [sweep_struct_master(s).sweep_info_neq, sweep_struct_master(s).sweep_results_neq] = ...
                    param_sweep_multi_v3([rate_index_sym sweep_index_vec(s)],...
                                        functionPath,sweep_options{:},...
                                        'half_max_flag',false,'equilibrium_flag',false,...
                                        'TauCycleTime',1,'downsample_output',1,'paramBounds',paramBounds); 

    [sweep_struct_master(s).sweep_info_eq, sweep_struct_master(s).sweep_results_eq] = ...
                    param_sweep_multi_v3([rate_index_sym sweep_index_vec(s)],...
                                        functionPath,sweep_options{:},...
                                        'half_max_flag',false,'equilibrium_flag',true,...
                                        'TauCycleTime',1,'downsample_output',1,'paramBounds',paramBounds); 
    toc
end                                  

%% Make scatter plots
close all
n_plot = 1500;
ylim_cell = {[0 sqrt(6e3)],[0 1.1],[0 0.55],[0 2.2]};
r_bins = linspace(0,1,50);
for s = 1:length(sweep_struct_master)
  
    % eq
    rate_vec_eq = sweep_struct_master(s).sweep_results_eq.metric_array(:,rate_index_sym);
    s_vec_eq = sweep_struct_master(s).sweep_results_eq.metric_array(:,sharpness_index_sym);
    metric_vec_eq = sweep_struct_master(s).sweep_results_eq.metric_array(:,sweep_index_vec(s));            
    options_eq = find(s_vec_eq>=0);
    [N_eq,~,bin_eq] = histcounts(rate_vec_eq(options_eq),r_bins);
    plot_indices_eq = randsample(options_eq,n_plot,true,1./N_eq(bin_eq));
    
    % neq
    rate_vec_neq = sweep_struct_master(s).sweep_results_neq.metric_array(:,rate_index_sym);
    s_vec_neq = sweep_struct_master(s).sweep_results_neq.metric_array(:,sharpness_index_sym);
    metric_vec_neq = sweep_struct_master(s).sweep_results_neq.metric_array(:,sweep_index_vec(s));
    options_neq = find(s_vec_neq>=0);
    [N_neq,~,bin_neq] = histcounts(rate_vec_neq(options_neq),r_bins);
    plot_indices_neq = randsample(options_neq,n_plot,true,1./N_neq(bin_neq));    
    
    if s <=2 
        metric_vec_eq = sqrt(exp(metric_vec_eq));
        metric_vec_neq = sqrt(exp(metric_vec_neq));
    end
    
    temp_fig = figure;
    cmap = brewermap(8,'Set2');
    

    % plot eq and neq domains
    hold on
    scatter(rate_vec_neq(plot_indices_neq),metric_vec_neq(plot_indices_neq),markerSize,'MarkerFaceColor',...
                cmap(2,:),'MarkerEdgeColor','k','MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5)
    scatter(rate_vec_eq(plot_indices_eq),metric_vec_eq(plot_indices_eq),markerSize,'MarkerFaceColor',...
                cmap(3,:),'MarkerEdgeColor','k','MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5)           
       

    xlabel('transcription rate (r)');
    ylabel(y_axis_cell{s})
    
    legend('non-equilibrium','equilibrium','Location','southeast')
    
    grid on
    set(gca,'FontSize',14)
%     set(gca,'Color',[228,221,209]/255) 

    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';

    
    temp_fig.InvertHardcopy = 'off';
    set(gcf,'color','w');

    box on
    
    ylim(ylim_cell{s})
%     xlim([0 2.2])

    saveas(temp_fig,[FigPath 'rate_vs_' sweep_id_cell{s} '.png'])
    saveas(temp_fig,[FigPath 'rate_vs_' sweep_id_cell{s} '.pdf'])
end
%% Plot illustrative titration curves

rate_bounds = [0.295 0.305];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find optimal neq systems
pd_vec_neq = sweep_struct_master(4).sweep_results_neq.metric_array(:,rate_index_sym);
sharp_vec_neq = sweep_struct_master(4).sweep_results_neq.metric_array(:,sharpness_index_sym);
sharp_raw_vec_neq = sweep_struct_master(4).sweep_results_neq.metric_array(:,sharpness_raw_index_sym);
prec_vec_neq = exp(sweep_struct_master(4).sweep_results_neq.metric_array(:,precision_index_sym));

pd_indices_neq = find(pd_vec_neq>=rate_bounds(1)&pd_vec_neq<=rate_bounds(2)&sharp_vec_neq>=1);

rate_array_sym_neq =- sweep_struct_master(4).sweep_results_neq.rate_array;

% overall optimum
[P,p_i_neq] = nanmax(prec_vec_neq(pd_indices_neq));
S = sharp_vec_neq(pd_indices_neq(p_i_neq));
s = sharp_raw_vec_neq(pd_indices_neq(p_i_neq));
s_rates_neq = rate_array_sym_neq(pd_indices_neq(p_i_neq),:);

% make arrays and calculate induction curves
c_vec = logspace(-2,2,1e3)';

s_rate_array_neq = repmat(s_rates_neq,length(c_vec),1);
s_rate_array_neq(:,1) = c_vec;

% map to correct subfolder
functionPath = sweep_struct_master(4).sweep_info_neq.functionPath;
addpath(genpath(functionPath));

% calculate production rates and noise profiles
paramCellNeq = mat2cell(s_rate_array_neq,size(s_rate_array_neq,1),ones(1,size(s_rate_array_neq,2)));

% get predicted response
r_curve_neq = productionRateFunction(paramCellNeq{:});

% now find corresponding Hill Function
Kd = 1*((1-0.3)/0.3)^(1/S);

r_curve_hill = c_vec.^S ./ (c_vec.^S + Kd^S);
%%
close all
S_fig = figure;

plot(c_vec,r_curve_neq,'-k','LineWidth',3)
set(gca,'xscale','log')
hold on
plot(c_vec,r_curve_hill,'--','Color',cmap(2,:),'LineWidth',3)
xlim([0.1 10])

% plot tangent curve
% tan_curve = s*c_vec;
% tan_curve = tan_curve-(mean(tan_curve(round(c_vec,2)==1))-0.3);
% plot(c_vec,tan_curve)
ylabel('transcription rate (r)');
xlabel('activator concentration (c)')

legend('input-ouptut function','equivalent Hill function','Location','northwest')

grid on
set(gca,'FontSize',14)

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

S_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

box on

ylim([0 1])
%     xlim([0 2.2])
saveas(S_fig,[FigPath 'example_curve.png'])
saveas(S_fig,[FigPath 'example_curve.pdf'])
%%

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