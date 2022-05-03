% script to generate illustrative log likelihood plots
clear
close all
addpath(genpath('../utilities'))

% make figure path
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\manuscript\Specificity\';
mkdir(DropboxFolder)

% Do parameter sweeps
% set basic parameters
nStates = 6;
rate_bounds = repmat([-6 ; 6],1,3*nStates-4); % constrain transition rate magnitude
[~,metric_names] = calculateMetricsSym([]);

% specify function path
functionPath = ['../utilities/metricFunctions/n' num2str(nStates) '_OR/'];

% make sure we're linked to the appropriate function subfolder% make sure we're linked to the appropriate function subfolder
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));

% get index of useful metrics
flux_index = find(strcmp(metric_names,'Flux'));
rate_index = find(strcmp(metric_names,'Production Rate'));
spec_index = find(strcmp(metric_names,'Specificity'));
spec_alt_index = find(strcmp(metric_names,'specFactorAlt'));
affinity_index = find(strcmp(metric_names,'AffinityVec'));

% set sim options
sweep_options = {'n_seeds',5,'n_iters_max',50,'nStates',nStates};

% sweep parallel specificity definition
tic
[sim_info_eq_alt, sim_struct_eq_alt] = param_sweep_multi_v2([affinity_index spec_alt_index], functionPath, sweep_options{:},...
                                          'half_max_flag',false,'cw',1,'equilibrium_flag',true);

[sim_info_eq, sim_struct_eq] = param_sweep_multi_v2([affinity_index spec_index], functionPath,sweep_options{:},...
                                          'half_max_flag',false,'cw',1,'equilibrium_flag',true);
                                        
[sim_info_neq, sim_struct_neq] = param_sweep_multi_v2([affinity_index spec_index], functionPath,sweep_options{:},...
                                          'half_max_flag',false,'cw',1,'equilibrium_flag',false);
toc
%% %%%%%%%%%%%%%%%%%%%%%%% define core parameters %%%%%%%%%%%%%%%%%%%%%%%%
n_plot = 3e3;
markerAlpha = 0.25; % marker transparency
markerSize = 75; % marker size
alpha = sim_info_eq_alt.specFactor;

% extract metric arrays
metric_array_eq = vertcat(sim_struct_eq.metric_array);
metric_array_eq_alt = vertcat(sim_struct_eq_alt.metric_array);
metric_array_neq = vertcat(sim_struct_neq.metric_array);

% pull affinity vectors (x axis)
affinity_vec_eq = 10.^metric_array_eq(:,affinity_index);
affinity_vec_eq_alt = 10.^metric_array_eq_alt(:,affinity_index);
affinity_vec_neq = 10.^metric_array_neq(:,affinity_index);

rate_vec_eq = metric_array_eq(:,rate_index);
rate_vec_eq_alt = metric_array_eq_alt(:,rate_index);
rate_vec_neq = metric_array_neq(:,rate_index);

% % pull specificity vectors
spec_vec_eq = 10.^metric_array_eq(:,spec_index);
spec_vec_eq_alt = 10.^metric_array_eq_alt(:,spec_alt_index);
spec_alt_filter = spec_vec_eq_alt>=0.01;
spec_vec_neq = 10.^metric_array_neq(:,spec_index);

% randomly select indices to plot
rng(475);
plot_indices_eq = randsample(1:length(spec_vec_eq),n_plot,false);
plot_indices_eq_alt = randsample(find(spec_alt_filter),n_plot,true);

% generate weights for neq
% neq_weights_1 = exp(-abs(spec_vec_neq - 2));
% neq_weights_2 = exp(-abs(spec_vec_neq - 0));

plot_indices_neq = randsample(1:length(spec_vec_neq),n_plot,false);

% make figure
close all
spec_comp_fig = figure;
cmap = brewermap(8,'Set2');
hold on

% % plot scatters
% scatter(affinity_vec_neq(plot_indices_neq),spec_vec_neq(plot_indices_neq),...
%           markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(2,:));     
% sneq = scatter(0,0,75,'MarkerEdgeColor','k', 'MarkerFaceColor',cmap(2,:));

scatter(affinity_vec_eq_alt(plot_indices_eq_alt),spec_vec_eq_alt(plot_indices_eq_alt),...
          markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(8,:));     
seq_alt = scatter(0,0,75,'MarkerEdgeColor','k', 'MarkerFaceColor',cmap(8,:));


% plot bounds
% plot(logspace(-3.5,3.5),repelem(2,50),'--','Color',cmap(2,:),'LineWidth',3)
ylim([-2.05 1])
xlim([10^-3.5 10^3.5])
 
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

precision_vs_sharpness.InvertHardcopy = 'off';
set(gcf,'color','w');

set(gca,'xscale','log')
set(gca,'yscale','log')
spec_comp_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
    
xlabel('affinity (k_{off}/k_{on})')
ylabel('activator fidelity (f/\alpha)')

ylim([10^-2 10^2.2])
xlim([10^-3.5 10^3.5])
set(gca,'ytick',[1e-2 1 1e2],'yticklabels',{'\alpha^{-1}','\alpha^{0}','\alpha^{1}'})

saveas(spec_comp_fig,[DropboxFolder 'affinity_par_spec.png'])
% saveas(spec_comp_fig,[DropboxFolder 'affinity_par_spec.pdf'])

scatter(affinity_vec_eq(plot_indices_eq),spec_vec_eq(plot_indices_eq),...
          markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(3,:));     
seq = scatter(0,0,75,'MarkerEdgeColor','k', 'MarkerFaceColor',cmap(3,:));

legend([seq seq_alt],'equilibrium (single locus)','equilibrium (parallel loci)','Location','northeast','Color','w')

saveas(spec_comp_fig,[DropboxFolder 'affinity_compare_spec.png'])
% saveas(spec_comp_fig,[DropboxFolder 'compare_spec.pdf'])

sneq1 = scatter(affinity_vec_neq(plot_indices_neq),spec_vec_neq(plot_indices_neq),...
          markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(2,:));     
sneq = scatter(0,0,75,'MarkerEdgeColor','k', 'MarkerFaceColor',cmap(2,:));

legend([seq_alt seq sneq],'parallel loci (eq.)','single locus (eq.)','single locus (non-eq.)','Location','southeast')
uistack(sneq1,'bottom')
saveas(spec_comp_fig,[DropboxFolder 'affinity_compare_spec_neq.png'])

%%
close all

spec_neq_fig = figure;
cmap = brewermap(8,'Set2');
hold on

% plot scatters
scatter(affinity_vec_neq(plot_indices_neq),spec_vec_neq(plot_indices_neq),...
          markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(2,:));     
sneq = scatter(0,0,75,'MarkerEdgeColor','k', 'MarkerFaceColor',cmap(2,:));

scatter(affinity_vec_eq(plot_indices_eq),spec_vec_eq(plot_indices_eq),...
          markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(3,:));     
seq = scatter(0,0,75,'MarkerEdgeColor','k', 'MarkerFaceColor',cmap(3,:));

scatter(affinity_vec_eq_alt(plot_indices_eq_alt),spec_vec_eq_alt(plot_indices_eq_alt),...
          markerSize*.75,'s','MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(8,:));     
seq_alt = scatter(0,0,75,'s','MarkerEdgeColor','k', 'MarkerFaceColor',cmap(8,:));

% plot bounds
% plot(logspace(-3.5,3.5),repelem(2,50),'--','Color',cmap(2,:),'LineWidth',3)
ylim([10^-2 10^2.2])
xlim([10^-3.5 10^3.5])

legend([seq_alt seq sneq],'parallel loci (eq.)','single locus (eq.)','single locus (non-eq.)','Location','southeast')

xlabel('sharpness (s)');
ylabel('precision (1/\sigma^2)')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

spec_neq_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

set(gca,'xscale','log')
set(gca,'yscale','log')

    
xlabel('effective K_D (k_{off}/k_{on})')
ylabel('activator fidelity (f/\alpha)')

saveas(spec_neq_fig,[DropboxFolder 'neq_spec.png'])
saveas(spec_neq_fig,[DropboxFolder 'neq_spec.pdf'])

%%
close all

spec_neq_fig = figure;
cmap = brewermap(8,'Set2');
hold on

% plot scatters
scatter(rate_vec_neq(plot_indices_neq),spec_vec_neq(plot_indices_neq),...
          markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(2,:));     
sneq = scatter(0,0,75,'MarkerEdgeColor','k', 'MarkerFaceColor',cmap(2,:));

scatter(rate_vec_eq(plot_indices_eq),spec_vec_eq(plot_indices_eq),...
          markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(3,:));     
seq = scatter(0,0,75,'MarkerEdgeColor','k', 'MarkerFaceColor',cmap(3,:));

scatter(rate_vec_eq_alt(plot_indices_eq_alt),spec_vec_eq_alt(plot_indices_eq_alt),...
          markerSize,'s','MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(8,:));     
seq_alt = scatter(0,0,75,'s','MarkerEdgeColor','k', 'MarkerFaceColor',cmap(8,:));

% plot bounds
% plot(logspace(-3.5,3.5),repelem(2,50),'--','Color',cmap(2,:),'LineWidth',3)
ylim([10^-2 10^2.2])
% xlim([10^-3.5 10^3.5])

legend([seq_alt seq sneq],'parallel loci (eq.)','single locus (eq.)','single locus (non-eq.)','Location','southwest')

xlabel('sharpness (s)');
ylabel('precision (1/\sigma^2)')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

spec_neq_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

set(gca,'yscale','log')

    
xlabel('fraction of time active (p_{on})')
ylabel('activator fidelity (f/\alpha)')

saveas(spec_neq_fig,[DropboxFolder 'neq_spec_pon.png'])
saveas(spec_neq_fig,[DropboxFolder 'neq_spec_pon.pdf'])

% % [affinity_vec, si] = sort(10.^sim_struct_eq(3).metric_array(:,affinity_index));

% 
% spec_parallel_vec = 10.^sim_struct_eq(3).metric_array(:,spec_alt_index);
% spec_parallel_vec = imgaussfilt(spec_parallel_vec(si),10);
% 
% spec_same_vec = 10.^sim_struct_eq(3).metric_array(:,spec_index);
% spec_same_vec = imgaussfilt(spec_same_vec(si),10);
% 
% plot_indices = sort(randsample(1:length(spec_same_vec),1e3,false));
% 
% [affinity_vec_neq, si_neq] = sort(10.^sim_struct_neq(3).metric_array(:,affinity_index));
% spec_same_vec_neq = 10.^sim_struct_neq(3).metric_array(:,spec_index);
% spec_same_vec_neq = spec_same_vec_neq(si_neq);
% 
% neq_aff_vec = logspace(-4,2,101); 
% aff_axis = 10.^(log10(neq_aff_vec(1:end-1)) + diff(log10(neq_aff_vec)));
% neq_spec_vec = NaN(size(aff_axis));
% 
% for i = 1:length(neq_aff_vec)-1
%     lb = neq_aff_vec(i);
%     ub = neq_aff_vec(i+1);
%     neq_spec_vec(i) = nanmax(spec_same_vec_neq(affinity_vec_neq>=lb & affinity_vec_neq < ub));
% end

% purple = [212 200 227]/255;
% 
% close all
% 
% specificity_fig = figure;
% cmap = brewermap([],'Set2');
% hold on
% 
% plot(affinity_vec(plot_indices),spec_parallel_vec(plot_indices),'Color','k','LineWidth',3)
% plot(affinity_vec(plot_indices),spec_same_vec(plot_indices),'Color',cmap(3,:),'LineWidth',3)
% plot(aff_axis,neq_spec_vec,'-.','Color',cmap(2,:),'LineWidth',3)
% 
% legend('f_p (equilibrium)','f_0 (equilibrium)','f_0 (non-eq bound)','Location','southeast')
% set(gca,'FontSize',14)
% %     set(gca, 'xtick', [],'ytick', [])
% xlabel('activator affinity (k_{off}/k_{on})')
% ylabel('activator fidelity (f_0/\beta)')
% set(gca,'Color',[228,221,209]/255) 
% ax = gca;
% 
% % ax.XAxis.MinorTickValues = 10.^(-5:.1:5);
% ax.YAxis(1).Color = 'k';
% ax.XAxis(1).Color = 'k';
% specificity_fig.InvertHardcopy = 'off';
% set(gcf,'color','w');
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% grid on
% 
% ylim([10^-2.5 10^2.5])
% xlim([1e-3 1e1])
% saveas(specificity_fig,[DropboxFolder 'compare_spec.png'])
% saveas(specificity_fig,[DropboxFolder 'compare_spec.pdf'])