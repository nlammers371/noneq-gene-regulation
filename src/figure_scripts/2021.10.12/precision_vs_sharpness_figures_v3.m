% Plot results precision vs. sharpness parameter sweeps
clear 
close all
addpath(genpath('../utilities/'))

nStates = 4;
functionPath = ['../utilities/metricFunctions/n' num2str(nStates) '_OR/'];

% make sure we're linked to the appropriate function subfolder
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));

% %%%%%%%%%%%%%%%%  set relative read and write paths %%%%%%%%%%%%%%%%%%%%
DataPath = ['../../out/bivariate_parameter_sweeps_n4_OR' filesep];
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\manuscript\';
FigPath = [DropboxFolder 'precision_vs_sharpness' filesep];
mkdir(FigPath);

% %%%%%%%%%%%%%%%%  get metric names %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,metric_names] = calculateMetricsSym([]);
flux_index = find(strcmp(metric_names,'Flux'));
rate_index = find(strcmp(metric_names,'Production Rate'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
precision_index = find(strcmp(metric_names,'Precision'));
information_rate_index = find(strcmp(metric_names,'DecisionRateNorm'));
decision_time_index = find(strcmp(metric_names,'DecisionTimeNorm'));
cycle_time_index = find(strcmp(metric_names,'CycleTime'));
phi_index = find(strcmp(metric_names,'Phi'));
state_entropy_index = find(strcmp(metric_names,'stateEntropy'));
rate_entropy_index = find(strcmp(metric_names,'rateEntropy'));

% %%%%%%%%%%%%%%%%  Set plot parameters and constants %%%%%%%%%%%%%%%%%%%%

n_plot = 3e3; % number of points to plot
markerAlpha = 0.5; % marker transparency
markerSize = 75; % marker size

% set sweep options
sweep_options = {'n_seeds',5,'n_iters_max',50,'nStates',nStates,'numCalcFlag',0};


%% %%%%%%%%%%%%%%%% Plot Sharpness vs Precision %%%%%%%%%%%%%%%%
% run sweep
tic
[sp_sim_info_neq, sim_struct_neq] = param_sweep_multi_v2([sharpness_index precision_index],functionPath,sweep_options{:},...
                                                            'half_max_flag',true,'equilibrium_flag',false,'saturationFlag',0);
                                                          
[sp_sim_info_eq, sp_sim_struct_eq] = param_sweep_multi_v2([sharpness_index precision_index],functionPath,sweep_options{:},...
                                                            'half_max_flag',true,'equilibrium_flag',true,'saturationFlag',0);
toc                                                          
                                                          
% process neq data
metric_array_neq = vertcat(sim_struct_neq.metric_array);   
sharpness_vec_neq = metric_array_neq(:,sharpness_index);
precision_vec_neq = exp(metric_array_neq(:,precision_index)).^2;
plot_filter_neq = sharpness_vec_neq>=0;% & precision_vec >=0;
plot_indices_neq = randsample(find(plot_filter_neq),n_plot,true);

% process eq data
metric_array_eq = vertcat(sp_sim_struct_eq.metric_array);   
sharpness_vec_eq = metric_array_eq(:,sharpness_index);
precision_vec_eq = exp(metric_array_eq(:,precision_index)).^2;
plot_filter_eq = sharpness_vec_eq>=0;% & precision_vec >=0;
plot_indices_eq = randsample(find(plot_filter_eq),n_plot,true);


cmap = brewermap(8,'Set2');
for i = 1:2
  
    precision_vs_sharpness = figure;
    
    hold on
    if i == 2
        scatter(sharpness_vec_neq(plot_indices_neq),precision_vec_neq(plot_indices_neq),...
              markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(2,:)); 
        sneq = scatter(0,0,75,'MarkerEdgeColor','k', 'MarkerFaceColor',cmap(2,:));
    end
    scatter(sharpness_vec_eq(plot_indices_eq),precision_vec_eq(plot_indices_eq),...
          markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(3,:));     
    seq = scatter(0,0,75,'MarkerEdgeColor','k', 'MarkerFaceColor',cmap(3,:));

    ylim([0 16.5])
    xlim([0 0.5])
        
    if i == 2
        legend([sneq seq],'nonequilibrium networks','equilibrium networks','Location','northeast')
    else
        legend([seq],'equilibrium networks','Location','northeast')
    end
    xlabel('sharpness (s)');
    ylabel('precision (1/\sigma^2)')
    grid on
    set(gca,'FontSize',14)
    set(gca,'Color',[228,221,209]/255) 

    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';

    precision_vs_sharpness.InvertHardcopy = 'off';
    set(gcf,'color','w');

    if i == 2          
        saveas(precision_vs_sharpness,[FigPath 'precision_vs_sharpness.png'])
        saveas(precision_vs_sharpness,[FigPath 'precision_vs_sharpness.pdf'])
    else
        saveas(precision_vs_sharpness,[FigPath 'precision_vs_sharpness_eq.png'])
        saveas(precision_vs_sharpness,[FigPath 'precision_vs_sharpness_eq.pdf'])
    end
end

%% Plot network heterogeneity measures
[sim_info_neq_end, sim_struct_neq_ent] = param_sweep_multi_v2([state_entropy_index rate_entropy_index],functionPath,sweep_options{:},...
                                                            'half_max_flag',true,'equilibrium_flag',false,'saturationFlag',0);

metric_array_ent = vertcat(sim_struct_neq_ent.metric_array);
rate_entropy_vec_full = metric_array_ent(:,rate_entropy_index);
state_entropy_vec_full = metric_array_ent(:,state_entropy_index);
sharpness_vec_ent = metric_array_ent(:,sharpness_index);

plot_filter_ent = sharpness_vec_ent>=0;% & precision_vec >=0;
sample_wts = (1./rate_entropy_vec_full(plot_filter_ent)) .* (1./state_entropy_vec_full(plot_filter_ent)).^.2;
sample_wts(sample_wts>10) = 10;
plot_indices_ent = randsample(find(plot_filter_ent),2*n_plot,true,sample_wts);

% extract metrics
sharpness_vec = metric_array_neq(:,sharpness_index);
precision_vec = metric_array_neq(:,precision_index);
info_vec = metric_array_neq(:,information_rate_index);

% extract rates and calculate state probabilities
rate_array_neq = vertcat(sim_struct_neq.rate_array);
rate_array_norm = rate_array_neq(:,2:end)./ sum(rate_array_neq(:,2:end),2);
paramCell = mat2cell(rate_array_neq,size(rate_array_neq,1),ones(1,size(rate_array_neq,2)));
stateProbArray = steadyStateVecFunction(paramCell{:});

% calculate rate and prob entropy measures
rate_entropy_vec = sum(rate_array_norm.*log(rate_array_norm),2);
rate_entropy_vec = rate_entropy_vec - min(rate_entropy_vec);
rate_entropy_vec = rate_entropy_vec/nanmax(rate_entropy_vec);

state_entropy_vec = sum(stateProbArray.*log(stateProbArray),2);
state_entropy_vec = state_entropy_vec - min(state_entropy_vec);
state_entropy_vec = state_entropy_vec/nanmax(state_entropy_vec);

% filter for high-performing networks
sharpness_vec(isnan(sharpness_vec))=-Inf;
[~,si_sharp] = sort(sharpness_vec,'descend');

precision_vec(isnan(precision_vec))=-Inf;
[~,si_precise] = sort(precision_vec,'descend');

% precise_val = prctile(precision_vec,99.9);
% info_val = prctile(info_vec,99.9);
% sharp_filter = sharpness_vec >= sharp_val;
% precise_filter = precision_vec >= precise_val;
% info_filter = info_vec >= info_val;
close all

topology_fig1 = figure;
cmap = brewermap(8,'Set2');
hold on

scatter(state_entropy_vec_full,rate_entropy_vec_full,...
          markerSize,'MarkerEdgeAlpha',.1,'MarkerEdgeColor','k','MarkerFaceAlpha',0.1, 'MarkerFaceColor',cmap(8,:));  
% scatter(state_entropy_vec(si_sharp(1:100)),rate_entropy_vec(si_sharp(1:100)),...
%           markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(4,:));     
% % ss = scatter(0,0,75,'MarkerEdgeColor','k', 'MarkerFaceColor',cmap(4,:));
% scatter(state_entropy_vec(si_precise(1:100)),rate_entropy_vec(si_precise(1:100)),...
%           markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(5,:));     

xlabel('state probability dispersion');
ylabel('transition rate dispersion')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
xlim([0 1])
ylim([0 1])
topology_fig1.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(topology_fig1,[FigPath 'dispersion_plot_bkg.png'])

topology_fig2 = figure;
cmap = brewermap(8,'Set2');
hold on

% scatter(state_entropy_vec_full,rate_entropy_vec_full,...
%           markerSize,'MarkerEdgeAlpha',.1,'MarkerEdgeColor','k','MarkerFaceAlpha',0.1, 'MarkerFaceColor',cmap(8,:));  
scatter(state_entropy_vec(si_sharp(1:100)),rate_entropy_vec(si_sharp(1:100)),...
          markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(4,:));     
ss = scatter(-1,-1,75,'MarkerEdgeColor','k', 'MarkerFaceColor',cmap(4,:));
scatter(state_entropy_vec(si_precise(1:100)),rate_entropy_vec(si_precise(1:100)),...
          markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(5,:));     
sp = scatter(-1,-1,75,'MarkerEdgeColor','k', 'MarkerFaceColor',cmap(5,:));

xlabel('state probability dispersion');
ylabel('transition rate dispersion')
legend([ss sp],'sharp networks','precise networks','Location','southeast')

grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
xlim([0 1])
ylim([0 1])

topology_fig2.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(topology_fig2,[FigPath 'dispersion_plot_main.png'])
saveas(topology_fig2,[FigPath 'dispersion_plot_main.pdf'])

%% Re-plot results from precision vs sharpness, this time interms of information gain
info_vec_neq = metric_array_neq(:,information_rate_index);
info_vec_eq = metric_array_eq(:,information_rate_index);

% generate normalized vectors
info_gain_vec = info_vec_neq ./ nanmax(info_vec_eq);
s_eq_max = nanmax(sharpness_vec_eq);
p_eq_max = nanmax(precision_vec_eq);
sharpness_gain_vec = (sharpness_vec_neq / s_eq_max).^2;
precision_gain_vec = precision_vec_neq / p_eq_max;

% calculate boundary for accessible region
nan_ft_neq = find(~isnan(sharpness_vec_neq)&sharpness_vec_neq>=0);
b_ids_neq = boundary(sharpness_vec_neq(nan_ft_neq),precision_vec_neq(nan_ft_neq),0.75);

sharp_bound_neq = sharpness_vec_neq(nan_ft_neq(b_ids_neq)).^2;
prec_bound_neq = precision_vec_neq(nan_ft_neq(b_ids_neq));

ds_vec = 1:5:length(prec_bound_neq);
prec_bound_neq = imgaussfilt(prec_bound_neq(ds_vec),2);
sharp_bound_neq = imgaussfilt(sharp_bound_neq(ds_vec),2);

% contour plot info 
i_max = 4;
levels = (0:(1/6):i_max)/(0.25^2*8);%max(info_array(:));
n_levels = floor(length(levels)/2);

% generate array of potential IR values
precision_vec = linspace(0,16);
sharpness_vec = linspace(0,0.5).^2;
contour_alpha = 0.6;
info_array = precision_vec'.*sharpness_vec/(s_eq_max^2*p_eq_max);

close all

info_gain_fig = figure;
cmap1 = flipud(brewermap(200,'Spectral'));
cmap2 = flipud(brewermap(66,'Spectral'));
cmap3 = vertcat(cmap2(1:50,:),cmap1(151:end,:));

colormap(cmap3)
hold on

% make contour plot first
[~,ch] = contourf(sharpness_vec/s_eq_max^2,precision_vec/p_eq_max,info_array,levels,'LineStyle','none');

% gray out inaccessible region
gray = [228,221,209]/255;
f1 = fillout(0,0,[0 4 0 2]);
f1.FaceColor = gray;
f1.FaceAlpha = 0.25;

f2 = fillout(sharp_bound_neq/s_eq_max^2,prec_bound_neq/p_eq_max,[0 4 0 2]);
f2.FaceColor = gray;%[0.5 0.5 0.5];
f2.FaceAlpha = 0.75;


% now overlay scatters
scatter(sharpness_gain_vec(plot_indices_neq),precision_gain_vec(plot_indices_neq), markerSize,info_gain_vec(plot_indices_neq),...
       'filled','MarkerEdgeAlpha',.9,'MarkerEdgeColor','k','MarkerFaceAlpha',0.9); 

h = colorbar; 
% caxis([0 4])
xlabel('sharpness gain (s^2/s^2_{eq})');
ylabel('precision gain (\sigma^2_{eq}/\sigma^2)')

ylabel(h,'non-equilibrium IR gain')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
h.Color = 'k';
info_gain_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
saveas(info_gain_fig,[FigPath 'info_gain_scatter.png'])
saveas(info_gain_fig,[FigPath 'info_gain_scatter.pdf'])

% %% Plot upper info bound, along with corresponding sharpness and precision
% n_points = 7;
% PhiMax = 15;
% phi_bins = [0 1e-3 1 2 5 15];
% 
% rng(123);
% tic
%                                                         
% [sim_info_neq, sim_struct_neq_ir] = param_sweep_multi_v2([phi_index information_rate_index],functionPath,sweep_options{:},...
%                                                             'half_max_flag',true,'equilibrium_flag',false);
% 
% toc   
%                                                 
% % track sharpness and precision of optimal networks Make figure
% metric_array_ir = vertcat(sim_struct_neq_ir.metric_array);
% ir_vec = metric_array_ir(:,information_rate_index);
% sharpness_vec = metric_array_ir(:,sharpness_index);
% precision_vec = metric_array_ir(:,precision_index);
% phi_vec = metric_array_ir(:,phi_index);
% 
% n_samples = 10;
% s_array = NaN(n_samples,length(phi_bins)-1);
% p_array = NaN(n_samples,length(phi_bins)-1);
% ir_array = NaN(n_samples,length(phi_bins)-1);
% 
% for m = 1:length(phi_bins)-1
%     
%     ir_temp = ir_vec;
%     ir_temp(isnan(sharpness_vec) | sharpness_vec < 0|phi_vec > phi_bins(m+1)) = -Inf;
%     [ir_sorted, mi_vec] = sort(ir_temp);
%     s_array(:,m) = sharpness_vec(mi_vec(end-n_samples+1:end)).^2;
%     p_array(:,m) = exp(precision_vec(mi_vec(end-n_samples+1:end))).^2;
% 
% end
% 
% % calculate eq and non-eq sharpness and precision boundaries
% nan_ft_neq = find(~isnan(sharpness_vec_neq)&sharpness_vec_neq>=0);
% b_ids_neq = boundary(sharpness_vec_neq(nan_ft_neq),precision_vec_neq(nan_ft_neq),0.5);
% 
% nan_ft_eq = find(~isnan(sharpness_vec_eq)&sharpness_vec_eq>=0);
% b_ids_eq = boundary(sharpness_vec_eq(nan_ft_eq),precision_vec_eq(nan_ft_eq),0.5);
% 
% % make figure
% 
% % downsample boundary contours
% plot_indices_neq = randsample(find(plot_filter_neq),3e4,false);
% % info_vec = sharpness_vec_neq.*precision_vec_neq;
% close all
% c_factor = ((sim_info_neq.cr1-sim_info_neq.cr0)/sim_info_neq.crs)^2;
% precision_vec = linspace(0,16);
% sharpness_vec = linspace(0,0.5).^2;
% contour_alpha = 0.6;
% info_array = precision_vec'.*sharpness_vec/(0.25^2*8);
% % n_levels = 7;
% i_max = 4;
% levels = (0:(1/6):i_max)/(0.25^2*8);%max(info_array(:));
% n_levels = floor(length(levels)/2);
% ir_axis = 0:1:i_max;
% ir_labels = round((0:.25:1) * c_factor * log2(exp(1)),3);
% % downsample 
% sharp_bound_eq = [0 0 0.0625 0.0625];
% prec_bound_eq = [0 8 8 0];
% s_eq = 0.25^2;
% p_eq = 8;
% 
% sharp_bound_neq = sharpness_vec_neq(nan_ft_neq(b_ids_neq)).^2;
% prec_bound_neq = precision_vec_neq(nan_ft_neq(b_ids_neq));
% 
% ds_vec = 1:5:length(prec_bound_neq);
% prec_bound_neq = imgaussfilt(prec_bound_neq(ds_vec),2);
% sharp_bound_neq = imgaussfilt(sharp_bound_neq(ds_vec),2);
% %%
% % close all
% contour_fig = figure;
% 
% hold on
% str1 = 'Spectral';
% str2 = 'Purples';
% cmap1 = (brighten(brewermap(length(phi_bins)-1,str1),0.25));
% cmap2 = contour_alpha*flipud(brewermap(n_levels+1,str2))+([1 1 1]-contour_alpha);
% cmap3 = contour_alpha*flipud(brewermap(n_levels^2,str2))+([1 1 1]-contour_alpha);
% % cmap2 = cmap2(5:end,:);
% cmap4 = [cmap2(1:end-1,:) ; cmap3(end-n_levels+1:end,:)];
% colormap(cmap4)
% [~,ch] = contourf(sharpness_vec/s_eq,precision_vec/p_eq,info_array,levels,'LineStyle','none');
% % imagesc(sharpness_vec/s_eq,precision_vec/p_eq,info_array)
% % gray out rest of plot region
% 
% f = fillout(sharp_bound_neq/s_eq,prec_bound_neq/p_eq,[0 4 0 2]);
% f.FaceColor = [228,221,209]/255;%[0.5 0.5 0.5];
% f.FaceAlpha = 0.5;
% 
% plot(sharp_bound_neq/s_eq,prec_bound_neq/p_eq,'-','Color','k','LineWidth',2)
% plot(sharp_bound_eq/s_eq,prec_bound_eq/p_eq,'--','Color','k','LineWidth',2)
% fill(sharp_bound_eq/s_eq,prec_bound_eq/p_eq,[228,221,209]/255,'FaceAlpha',0.2,'EdgeAlpha',0);
% 
% plot(mean(s_array)/s_eq, mean(p_array)/p_eq,'-','Color','k','LineWidth',3)
% for m = 1:length(phi_bins)-1
%     scatter(mean(s_array(:,m)/s_eq), mean(p_array(:,m)/p_eq),45,'MarkerFaceColor',cmap1(m,:),'MarkerEdgeColor','k',...
%               'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1)            
% end
% 
% 
% 
% xlabel('sharpness gain (s^2/s^2_{eq})');
% ylabel('precision gain (\sigma^2_{eq}/\sigma^2)')
% h = colorbar;
% h.Color = 'k';
% colormap(h,cmap2)
% % h.Ticks = 0:1:4;
% h.Ticks = 0:2:8;
% h.TickLabels = ir_labels;
% 
% % h.TickLabels = round(logspace(-2, log10(15), 5),2);
% % ylabel(h,'\Phi (k_bT per cycle)')
% ylabel(h,'IR (bits per cycle)')
% grid on
% set(gca,'FontSize',14)
% set(gca,'Color',[228,221,209]/255) 
% 
% ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.XAxis(1).Color = 'k';
% 
% saveas(contour_fig,[FigPath 'info_trajectory_fig.png'])
% saveas(contour_fig,[FigPath 'info_trajectory_fig.pdf'])
% 
% %% generate predicted response profiles for precision and sharpness motifs
% n_cycles = 100;
% 
% % [sp_sim_info_neq, sim_struct_neq] = param_sweep_multi_v2([sharpness_index precision_index],functionPath,sweep_options{:},...
% %                                                             'half_max_flag',true,'equilibrium_flag',false,'saturationFlag',1);
% 
% metric_array_neq = vertcat(sim_struct_neq.metric_array);
% rate_array_neq = vertcat(sim_struct_neq.rate_array);
% sharpness_vec = metric_array_neq(:,sharpness_index);
% precision_vec = metric_array_neq(:,precision_index);
% 
% [~,sharpness_i] = nanmax(sharpness_vec.*(exp(precision_vec).^2<=1));
% [~,precision_i] = nanmax(precision_vec.*(sharpness_vec>=.13));
% 
% sharp_parameters = rate_array_neq(sharpness_i,:);
% precise_parameters = rate_array_neq(precision_i,:);
% 
% % generate predicted ressponse profiles
% c_vec = logspace(-3,3,76);
% 
% % generate input cell arrays
% sharpMat = [c_vec' repmat(sharp_parameters(:,2:end),length(c_vec),1)];
% sharpParamCell = mat2cell(sharpMat,size(sharpMat,1),ones(1,size(sharpMat,2)));
% 
% preciseMat = [c_vec' repmat(precise_parameters(:,2:end),length(c_vec),1)];
% preciseParamCell = mat2cell(preciseMat,size(preciseMat,1),ones(1,size(preciseMat,2)));
% 
% % calculate production rate         
% sharp_response = productionRateFunction(sharpParamCell{:});
% sharp_var = intrinsicVarianceFunction(sharpParamCell{:});
% sharp_tOn = TauONFunction(sharpParamCell{:});
% sharp_tOff = TauOFFFunction(sharpParamCell{:});
% sharp_tau = sharp_tOn + sharp_tOff;
% sharp_noise = sqrt(sharp_var./sharp_tau / n_cycles);
%     
% s_ub = sharp_response + sharp_noise;
% % s_ub(s_ub>1) = 1;
% s_lb = sharp_response - sharp_noise;
% 
% precise_response = productionRateFunction(preciseParamCell{:});
% precise_var = intrinsicVarianceFunction(preciseParamCell{:});
% precise_tOn = TauONFunction(preciseParamCell{:});
% precise_tOff = TauOFFFunction(preciseParamCell{:});
% precise_tau = precise_tOn + precise_tOff;
% precise_noise = sqrt(precise_var./precise_tau / n_cycles);
% 
% p_ub = precise_response + precise_noise;
% p_lb = precise_response - precise_noise;
% 
% % make figure
% close all
% profile_fig = figure;
% hold on
% 
% fill([c_vec fliplr(c_vec)],[s_ub',fliplr(s_lb')],cmap(4,:),'FaceAlpha',0.2)
% fill([c_vec fliplr(c_vec)],[p_ub',fliplr(p_lb')],cmap(5,:),'FaceAlpha',0.2)
% p1 = plot(c_vec,sharp_response,'Color',cmap(4,:),'LineWidth',1.5);
% p2 = plot(c_vec,precise_response,'Color',cmap(5,:),'LineWidth',1.5);
% grid on
% set(gca,'xScale','log');
% 
% xlabel('activator concentration (c)');
% ylabel('transcription rate (r)')
% legend([p1 p2],'sharp gene circuit','precise gene circuit','Location','northwest')
% grid on
% set(gca,'FontSize',14)
% set(gca,'Color',[228,221,209]/255) 
% 
% ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.XAxis(1).Color = 'k';
% xlim([1e-2 1e2]);
% ylim([0 1])
% profile_fig.InvertHardcopy = 'off';
% set(gcf,'color','w');
% 
% saveas(profile_fig,[FigPath 'motif_response_profiles.png'])
% saveas(profile_fig,[FigPath 'motif_response_profiles.pdf'])
