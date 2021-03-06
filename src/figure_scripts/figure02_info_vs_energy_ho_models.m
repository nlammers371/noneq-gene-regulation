% Plot results for IR vs energy for higher order models
clear 
close all
addpath(genpath('../utilities/'))

% %%%%%%%%%%%%%%%%  set relative read and write paths %%%%%%%%%%%%%%%%%%%%
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\';
DataPath = [DropboxFolder  'SweepOutput\sweeps01_info_vs_energy\' filesep ];
FigPath = [DropboxFolder '\manuscript\info_vs_energy' filesep];
mkdir(FigPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%% Examine impact of adding binding sites %%%%%%%%%%%%%%%%%

% get metric names for numeric sweeps
[~,~,metric_names_num] = calculateMetricsNumeric_v3([]);

ir_index_num = find(strcmp(metric_names_num,'IR'));
rate_index_num = find(strcmp(metric_names_num,'ProductionRate'));
sharp_index_num = find(strcmp(metric_names_num,'Sharpness'));
phi_index_num = find(strcmp(metric_names_num,'Phi'));
inv_dtime_index_num = find(strcmp(metric_names_num,'InverseDecisionTime'));

% get list of sweep results files with only 1 genera TF reaction
multi_bs_sweep_files = dir([DataPath 'sweep_results*g01*']);
multi_bs_info_files = dir([DataPath 'sweep_info*g01*']);

% load
master_struct_multi_bs = struct;
for f = 1:length(multi_bs_sweep_files)
  
    load([DataPath multi_bs_sweep_files(f).name])
    load([DataPath multi_bs_info_files(f).name])
    
    master_struct_multi_bs(f).sweep_results = sim_results;
    master_struct_multi_bs(f).sweep_info = sim_info;
end    


%% plot Phi vs IR 

% calculate upper IR bounds for each NBS
phi_bins = logspace(-2.34,log10(75),201);
phi_axis = phi_bins(2:end);
ir_bound_array = NaN(length(phi_axis),length(master_struct_multi_bs));
ir_rates_array = NaN(length(phi_axis),length(master_struct_multi_bs),size(master_struct_multi_bs(end).sweep_results.rate_array,2));

for i = length(master_struct_multi_bs):-1:1
    metric_array = master_struct_multi_bs(i).sweep_results.metric_array;
    rate_array = master_struct_multi_bs(i).sweep_results.rate_array;
    phi_vec = metric_array(:,phi_index_num);
    ir_vec = metric_array(:,ir_index_num);
    r_vec = metric_array(:,rate_index_num);
    s_vec = metric_array(:,sharp_index_num);
    
    for p = 1:length(phi_bins)-1
      
        phi_filter = phi_vec>phi_bins(p) & phi_vec<=phi_bins(p+1);
        ir_filtered = ir_vec(phi_filter);
        rates_filtered = rate_array(phi_filter,:);
        
        % deal with outliers due to numerical precision issues (NL: need to
        % resolve this)
%         if phi_bins(p+1) <= 0.1
%            ir_filtered(ir_filtered > 0.067) = NaN;
%         end
        if any(ir_filtered)
            [ir_bound_array(p,i), mi] = nanmax(ir_filtered*log2(exp(1))); % convert to bits
            ir_rates = rates_filtered(mi,:);
            ir_rates_array(p,i,1:length(ir_rates)) = ir_rates;
        end
    end
end

close all
ds_factor = 3;
color_ind = 5;
sym_list = {'o','^','d','s','v'};
% bs_vec = 1:length(master_struct_multi_bs);

% set plot parameters
markerSize = 50;
rng(231);
alphaFactor = 0.25;
 
% make figure

% Define colormaps for use throughout
cmap_pu = brewermap(8,'Purples');
cmap_rd = brewermap(8,'Reds');
cmap_bu = brewermap(8,'Blues');
cmap_gre = brewermap(8,'Greens');
cmap_gra = brewermap(8,'Greys');

close all

cmap = [cmap_pu(3,:); cmap_gre(color_ind,:); cmap_rd(color_ind,:); ...
          cmap_bu(color_ind,:) ; cmap_gra(color_ind,:)];

         
% generate down-sampled version of x axis for plotting
ir_bound_array_fit = NaN(length(phi_axis),size(ir_bound_array,2));

% fit smoothing splines to upper bounds for each genecircuit type

for i = length(master_struct_multi_bs):-1:1  
    % fit smoothing spline
    ir_bound_temp = ir_bound_array(:,i);
    ir_bound_temp = ir_bound_temp(~isnan(ir_bound_temp));
    phi_temp = phi_axis(~isnan(ir_bound_array(:,i)));
    phi_max = max(phi_temp);
    phi_bs_axis_temp = phi_axis(phi_axis<=phi_max);
    
    [f_spline,~,~] = fit(log10(phi_temp'),ir_bound_temp,'smoothingspline','SmoothingParam',0.99);
    f_spline_ds = feval(f_spline,log10(phi_bs_axis_temp));%downsample(f_spline,ds_factor);    

    % plot smoothed trends
    ir_bound_array_fit(phi_axis<=phi_max,i) = f_spline_ds; 
end    

% remove phi values for which one more more IR values are missing
nan_flags = any(isnan(ir_bound_array_fit),2);
ir_bound_array_fit = ir_bound_array_fit(~nan_flags,:);
phi_axis = phi_axis(~nan_flags);

% s = [];

info_bs_fig = figure;

hold on
% plot area vectors
for i = length(master_struct_multi_bs):-1:1    
                      
    if i > 1
        fill([phi_axis fliplr(phi_axis)],[ir_bound_array_fit(:,i)' fliplr(ir_bound_array_fit(:,i-1)')],...
                                        cmap(i,:),'FaceAlpha',alphaFactor*3,'EdgeAlpha',0);
        plot(phi_axis,ir_bound_array_fit(:,i)','Color', brighten(cmap(i,:),-0.5),'LineWidth',3);
    else
       fill([phi_axis fliplr(phi_axis)],[ir_bound_array_fit(:,i)' repelem(1e-6,length(ir_bound_array_fit(:,i)))],...
                                        cmap(i,:),'FaceAlpha',alphaFactor*3,'EdgeAlpha',0);
       plot(phi_axis,ir_bound_array_fit(:,i)','Color', brighten(cmap(i,:),-0.5),'LineWidth',3);                                      
    end
end

% % overlay line plots
% for i = 1:length(bs_vec)
%     s(i) = scatter(phi_axis,ir_bound_array_fit(:,i),markerSize','MarkerFaceColor',brighten(cmap(i,:),-0.25),'MarkerEdgeColor','k');
% end    

ylim([1e-3 0.5])
xlim([5e-3 3e1])

% h = colorbar;
% ylabel(h,'number of binding sites')
% legend(s,'1 bs','2 bs','3 bs','4 bs','5 bs','Location','southeast','Orientation','horizontal');
xlabel('energy dissipation rate (k_BT per burst cycle)');
ylabel('information rate (bits per burst cycle)')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
set(gca,'xtick',[0.01 0.1 1 10])
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

info_bs_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
set(gca,'xscale','log')
set(gca,'yscale','log')

saveas(info_bs_fig,[FigPath 'info_vs_phi_bs.png'])
saveas(info_bs_fig,[FigPath 'info_vs_phi_bs.pdf']) 

%% Plot IR eq bounds as a function of N_b
% note that I'm using low Phi as a proxy for eq performance

n_plot = 500;

eq_bound_bs = figure;

hold on


% make plots
for i = 1:length(master_struct_multi_bs)
    ir_vec = master_struct_multi_bs(i).sweep_results.metric_array(:,ir_index_num);
    phi_vec = master_struct_multi_bs(i).sweep_results.metric_array(:,phi_index_num);
    num_err_bound = 1;
%     if i == 4
%         num_err_bound = ir_bound_array_fit(1,4)/log2(exp(1));
%     elseif i == 5
%         num_err_bound = 0.067;
%     end
    ir_options = find(phi_vec<=phi_bins(2) & ir_vec >=0 & ir_vec < num_err_bound);
    [N,~,bin] = histcounts(ir_vec(ir_options));
    plot_ids = randsample(ir_options,n_plot,true,1./N(bin));
    g_temp = i + normrnd(0,0.025,n_plot,1);
    scatter(g_temp,ir_vec(plot_ids)*log2(exp(1)),25,'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.1);
end

% plot quadratic trend
bs_vec = 0:6;
trend_vec = bs_vec.^2*ir_bound_array_fit(1,1);
plot(bs_vec,trend_vec,'--k','LineWidth',3);


% ylim([1e-3 0.5])
% legend([seq sneq],'equilibrium','non-equilibrium','Location','northwest')
xlabel('number of binding sites (N_b)');
ylabel('information rate (bits per burst cycle)')

grid on
box on

set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
ylim([0 0.1])
xlim([0.5 5.5])
eq_bound_bs.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(eq_bound_bs,[FigPath 'info_vs_bs.png'])
saveas(eq_bound_bs,[FigPath 'info_vs_bs.pdf']) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot DT eq bounds as a function of N_b
% note that I'm using low Phi as a proxy for eq performance

n_plot = 500;
eq_bound_bs = figure;

hold on

% make plots
for i = 1:length(master_struct_multi_bs)
    ir_vec = master_struct_multi_bs(i).sweep_results.metric_array(:,ir_index_num);
    phi_vec = master_struct_multi_bs(i).sweep_results.metric_array(:,phi_index_num);
    dt_vec = 1./master_struct_multi_bs(i).sweep_results.metric_array(:,inv_dtime_index_num );
    
    num_err_bound = 1;
%     if i == 4
%         num_err_bound = ir_bound_array_fit(1,4)/log2(exp(1));
%     elseif i == 5
%         num_err_bound = 0.067;
%     end
    dt_options = find(phi_vec<=phi_bins(2) & ir_vec >=0 & ir_vec < num_err_bound & ~isnan(dt_vec) & dt_vec <999 & dt_vec >= 0);
    [N,~,bin] = histcounts(dt_vec(dt_options),logspace(-1,3));
    plot_ids = randsample(dt_options,n_plot,true,1./N(bin));
    g_temp = i + normrnd(0,0.025,n_plot,1);
    scatter(g_temp,dt_vec(plot_ids),25,'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.1);
    if i == 1
        dt_min = nanmin(dt_vec(dt_options));
    end
end

% plot quadratic trend
bs_vec = linspace(0,6);
trend_vec = bs_vec.^-2*dt_min;
plot(bs_vec,trend_vec,'-.k','LineWidth',3);


% ylim([1e-3 0.5])
% legend([seq sneq],'equilibrium','non-equilibrium','Location','northwest')
xlabel('number of binding sites');
ylabel('decision time (burst cycles)')

grid on
box on

set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
ylim([1 5e2])
xlim([0.5 5.5])
set(gca,'yscale','log')
eq_bound_bs.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(eq_bound_bs,[FigPath 'dt_vs_bs_eq.png'])
saveas(eq_bound_bs,[FigPath 'dt_vs_bs_eq.pdf']) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%% Look at impact of multiple Activation Steps %%%%%%%%%%%%%%%

% get list of sweep results files with only 1 genera TF reaction
multi_g_sweep_files = dir([DataPath 'sweep_results_s01_ns00_g0*']);
multi_g_info_files = dir([DataPath 'sweep_info_s01_ns00_g0*']);
  
% load
master_struct_multi_g = struct;
for f = 1:length(multi_g_sweep_files)
  
    load([DataPath multi_g_sweep_files(f).name])
    load([DataPath multi_g_info_files(f).name])
    
    master_struct_multi_g(f).sweep_results = sim_results;
    master_struct_multi_g(f).sweep_info = sim_info;
end    

% calculate upper IR bounds

% calculate upper IR bounds for each NBS
phi_bins_g = logspace(-2.5,log10(1e4),201);
phi_axis_g = phi_bins_g(2:end);
ir_gen_array = NaN(length(phi_axis_g),length(master_struct_multi_g));
ir_gen_rates_array = NaN(length(phi_axis_g),length(master_struct_multi_g),size(master_struct_multi_g(end).sweep_results.rate_array,2));

for i = length(master_struct_multi_g):-1:1
    metric_array = master_struct_multi_g(i).sweep_results.metric_array;
    rate_array = master_struct_multi_g(i).sweep_results.rate_array;
    phi_vec = metric_array(:,phi_index_num);
    ir_vec = metric_array(:,ir_index_num);
    sharp_vec = metric_array(:,sharp_index_num);
    for p = 1:length(phi_bins_g)-1
        phi_filter = phi_vec>phi_bins_g(p) & phi_vec<=phi_bins_g(p+1);
        ir_filtered = ir_vec(phi_filter);
        rates_filtered = rate_array(phi_filter,:);
        if any(ir_filtered)            
            [ir_gen_array(p,i), mi] = nanmax(ir_filtered*log2(exp(1))); % convert to bits
            ir_rates = rates_filtered(mi,:);
            ir_gen_rates_array(p,i,1:length(ir_rates)) = ir_rates;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make figure 

% downsample phi axis
phi_axis_ds_g = downsample(phi_axis_g,ds_factor);
phi_axis_ds_g = unique([phi_axis_ds_g max(phi_axis_g)]);
ir_gen_array_sm = NaN(length(phi_axis_ds_g),size(ir_gen_array,2));

% fit upper bounds to each model
for i = length(master_struct_multi_g):-1:1  
  
    % fit spline to generate smoothed trend
    ir_bound_raw = ir_gen_array(:,i);            
    nan_filter = ~isnan(ir_bound_raw);
    ir_bound_temp = ir_bound_raw(nan_filter);
    phi_temp = phi_axis_g(nan_filter);
    phi_max = max(phi_temp);
    phi_axis_ds_temp = phi_axis_ds_g(phi_axis_ds_g<=phi_max);
    
    [f_spline,~,~] = fit(log10(phi_temp'),ir_bound_temp,'smoothingspline','SmoothingParam',0.99);
    f_spline_ds = feval(f_spline,log10(phi_axis_ds_temp));
    
    % store smoothed trends
    ir_gen_array_sm(phi_axis_ds_g<=phi_max,i) = f_spline_ds;
end    
 
close all

info_g_fig = figure;
s = [];
% colormap(cmap(4:end,:));
hold on

for i = size(ir_gen_array_sm,2):-1:1
                      
    ir_bound = ir_gen_array_sm(:,i);
    nan_filter = ~isnan(ir_bound);
    phi_axis = phi_axis_ds_g(nan_filter);
    ir_plot = ir_bound(nan_filter);
    
    if i > 2
        ir_bound_prev = ir_gen_array_sm(:,i-1);
        ir_prev = ir_bound_prev(nan_filter);
        fill([phi_axis fliplr(phi_axis)],[ir_plot' fliplr(ir_prev')],...
                                        cmap_pu(2+i,:),'FaceAlpha',0.9,'EdgeAlpha',0);
        plot(phi_axis,ir_plot','Color', brighten(cmap_pu(2+i,:),-0.25),'LineWidth',3);
    else
       fill([phi_axis fliplr(phi_axis)],[ir_plot' repelem(1e-6,length(ir_plot))],...
                                        cmap_pu(2+i,:),'FaceAlpha',0.9,'EdgeAlpha',0);
       plot(phi_axis,ir_plot','Color', brighten(cmap_pu(2+i,:),-0.25),'LineWidth',3);                                      
    end
end

% generate predicted bound
% ir_eq = nanmax(ir_gen_array(phi_axis_g<5e-3,1));
% ir_bound_pd = ir_eq + ir_eq*sqrt(phi_axis_g);
% plot(phi_axis_g,ir_bound_pd,'-k','LineWidth',4);  

set(gca,'xscale','log')
ylim([0 0.06])
xlim([5e-3 5e3])

% h = colorbar;
% ylabel(h,'number of binding sites')
% legend(s,'1 general TF','2 general TFs','3 general TFs','Location','northwest');
xlabel('energy dissipation rate (k_BT per burst cycle)');
ylabel('information rate (bits per burst cycle)')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
set(gca,'xtick',[1e-2 1e-1 1 1e1 1e2 1e3])
info_g_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
set(gca,'xscale','log')
% set(gca,'yscale','log')

saveas(info_g_fig,[FigPath 'info_vs_phi_g.png'])
saveas(info_g_fig,[FigPath 'info_vs_phi_g.pdf']) 


%% Plot IR neq bounds as a function of N_LC
% note that I'm using low Phi as a proxy for eq performance

n_plot = 500;

eq_bound_a = figure;
hold on
% make plots
ir_maxima = [];
for i = 1:length(master_struct_multi_g)
    ir_vec = master_struct_multi_g(i).sweep_results.metric_array(:,ir_index_num);
    phi_vec = master_struct_multi_g(i).sweep_results.metric_array(:,phi_index_num);
    
    
    ir_maxima(i) = nanmax(ir_vec)*log2(exp(1));
    
    ir_options = find(ir_vec >=0 & ir_vec);    
    [N,~,bin] = histcounts(ir_vec(ir_options));
    plot_ids = randsample(ir_options,n_plot,true,1./N(bin));
    g_temp = 1 + i + normrnd(0,0.025,n_plot,1);
    scatter(g_temp,ir_vec(plot_ids)*log2(exp(1)),25,'MarkerFaceColor',cmap_pu(2+i,:)...
                 ,'MarkerEdgeColor','k','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.1);
end

% plot quadratic trend
a_vec = 0:6;
mdl = fitlm(1+(0:4)',[0; ir_maxima']);%,'Intercept',false);
trend_vec = predict(mdl,a_vec');
plot(a_vec,trend_vec,'-.k','LineWidth',3);


% ylim([1e-3 0.5])
% legend([seq sneq],'equilibrium','non-equilibrium','Location','northwest')
xlabel('number of locus conformation (N_{LC})');
ylabel('information rate (bits per burst cycle)')

grid on
box on

set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
ylim([0 0.1])
xlim([1.5 5.5])
eq_bound_a.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(eq_bound_a,[FigPath 'info_vs_a.png'])
saveas(eq_bound_a,[FigPath 'info_vs_a.pdf']) 


%% Plot DT neq bounds as a function of N_LC
% note that I'm using low Phi as a proxy for eq performance

n_plot = 500;

neq_bound_a = figure;

hold on
dt_minima = [];
% make plots
for i = 1:length(master_struct_multi_g)
    ir_vec = master_struct_multi_g(i).sweep_results.metric_array(:,ir_index_num);
    phi_vec = master_struct_multi_g(i).sweep_results.metric_array(:,phi_index_num);
    dt_vec = 1./master_struct_multi_g(i).sweep_results.metric_array(:,inv_dtime_index_num );
       
    dt_options = find(ir_vec >=0 & ~isnan(dt_vec) & dt_vec <999 & dt_vec >= 0);
    [N,~,bin] = histcounts(dt_vec(dt_options),logspace(-1,3));
    plot_ids = randsample(dt_options,n_plot,true,1./N(bin));
    g_temp = i + 1 + normrnd(0,0.025,n_plot,1);
    
    scatter(g_temp,dt_vec(plot_ids),25,'MarkerFaceColor',cmap_pu(2+i,:),'MarkerEdgeColor',...
                                    'k','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.1);        
                                  
    dt_minima(i) = nanmin(dt_vec(dt_options));
    
end

% plot quadratic trend
a_vec = linspace(0,6);
mdl = fitlm(1+(1:4)',1./dt_minima');
trend_vec = predict(mdl,a_vec');
plot(a_vec,1./trend_vec,'-.k','LineWidth',3);

xlabel('number of locu conformation (N_{LC})');
ylabel('decision time (burst cycles)')

grid on
box on

set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
ylim([1 5e2])
xlim([1.5 5.5])
set(gca,'yscale','log')
neq_bound_a.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(neq_bound_a,[FigPath 'dt_vs_a_neq.png'])
saveas(neq_bound_a,[FigPath 'dt_vs_a_neq.pdf'])

% %%
% phi_min = nanmin(phi_vec(dt_vec<=7.9))

% %% Plot maximum IR rate vs Ng
% ng_ir_vec_neq = nanmax(ir_gen_array);
% ng_ir_vec_eq = nanmax(ir_gen_array(phi_axis_g<=5e-3,:));
% ng_vec = 1:4;
% % generate predicted neq trend
% ng_vec_interp = linspace(0,4);
% mdl_ng = fitlm(ng_vec,ng_ir_vec_neq);
% ng_ir_pd_neq = predict(mdl_ng,ng_vec_interp')';
% ng_ir_pd_eq = repelem(ng_ir_vec_eq(1),length(ng_vec_interp));
% 
% eq_neq_bound_ng = figure;
% hold on
% % plot eq and noneq regimes
% bottomline = repelem(-.1,length(ng_vec_interp));
% fill([ng_vec_interp fliplr(ng_vec_interp)], [ng_ir_pd_neq fliplr(ng_ir_pd_eq)],cmap2(2,:),'FaceAlpha',0.25,'EdgeAlpha',0)
% fill([ng_vec_interp fliplr(ng_vec_interp)], [ng_ir_pd_eq bottomline],cmap2(3,:),'FaceAlpha',0.5,'EdgeAlpha',0)
% 
% % make plots
% plot(ng_vec_interp,ng_ir_pd_eq,'-','Color','k','LineWidth',2)
% plot(ng_vec_interp,ng_ir_pd_neq,'-','Color','k','LineWidth',2)
% 
% % make scatters
% seq = scatter(ng_vec,ng_ir_vec_eq,75,'MarkerFaceColor',cmap2(3,:),'MarkerEdgeColor','k');
% sneq = scatter(ng_vec,ng_ir_vec_neq,75,'MarkerFaceColor',cmap2(2,:),'MarkerEdgeColor','k');
% 
% 
% legend([seq sneq],'equilibrium','non-equilibrium','Location','northwest')
% xlabel('number of general TF reactions');
% ylabel('information rate (bits per burst cycle)')
% 
% grid on
% box on
% 
% set(gca,'FontSize',14)
% set(gca,'Color',[228,221,209]/255) 
% 
% ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.XAxis(1).Color = 'k';
% ylim([0 0.07])
% xlim([1 4])
% eq_neq_bound_ng.InvertHardcopy = 'off';
% set(gcf,'color','w');
% 
% saveas(eq_neq_bound_ng,[FigPath 'info_vs_ng.png'])
% saveas(eq_neq_bound_ng,[FigPath 'info_vs_ng.pdf']) 
% 
% %% Do the same for decision times
% % initialize arrays to store rates
% % rate_array = NaN(n_lines,size(rate_array_flux,2));
% eq_ind_full_ng = [];
% neq_ind_full_ng = [];
% 
% for i = 1:length(ng_vec)
%     [~,eq_ind_full_ng(i)] = nanmax(ir_gen_array(phi_axis<=5e-3&phi_axis<=3e-3,i));
%     [~,neq_ind_full_ng(i)] = nanmax(ir_gen_array(:,i));
% end    
% 
% tau_array_neq_ng = NaN(1,length(ng_vec));
% tau_array_eq_ng = NaN(1,length(ng_vec));
% 
% 
% % iterate
% for i = 1:length(ng_vec)
%        
%   optimal_rates_neq = reshape(ir_gen_rates_array(neq_ind_full_ng(i),i,:),1,[]);
%   optimal_rates_neq = optimal_rates_neq(~isnan(optimal_rates_neq));
%   optimal_rates_eq = reshape(ir_gen_rates_array(eq_ind_full_ng(i),i,:),1,[]);
%   optimal_rates_eq = optimal_rates_eq(~isnan(optimal_rates_eq));
% 
%   % set function path
%   sweep_info_temp = master_struct_multi_g(i).sweep_info;
%   functionPath = sweep_info_temp.functionPath;
%   slashes = strfind(functionPath,'\');
%   simName = functionPath(slashes(end-1)+1:slashes(end)-1);
%   rmpath(genpath('../utilities/metricFunctions/'));
%   addpath(genpath(['../utilities/metricFunctions/numeric/' simName]));
% 
%   sweep_info_temp.minError = minError;
%   sweep_info_temp.numerical_precision = 5;
%   
%   tau_neq = NaN;
%   iter = 1;
%   while isnan(tau_neq) && iter < n_tries
%       metric_vec_neq = calculateMetricsNumeric_v3(optimal_rates_neq, sweep_info_temp);
%       tau_neq = 1./metric_vec_neq(inv_dtime_index_num);
%       iter = iter + 1;
%   end
%   tau_array_neq_ng(i) = tau_neq;
%   
%   tau_eq = NaN;
%   iter = 1;
%   while isnan(tau_eq) && iter < n_tries
%       metric_vec_eq = calculateMetricsNumeric_v3(optimal_rates_eq, sweep_info_temp);
%       tau_eq = 1./metric_vec_eq(inv_dtime_index_num);
%       iter = iter + 1;
%   end
%   tau_array_eq_ng(i) = tau_eq;
%   
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Make decision time vs binding site figure
% 
% topline = repelem(1e5,length(ng_vec));
% 
% eq_neq_dtime_ng = figure;
% cmap2 = brewermap(8,'Set2');
% 
% hold on
% 
% % plot eq and noneq regimes
% fill([ng_vec fliplr(ng_vec)], [tau_array_eq_ng topline],cmap2(3,:),'FaceAlpha',0.5,'EdgeAlpha',0)
% fill([ng_vec fliplr(ng_vec)], [tau_array_neq_ng fliplr(tau_array_eq_ng)],cmap2(2,:),'FaceAlpha',0.5,'EdgeAlpha',0)
% 
% % make plots
% % eq_tau_pd = tau_array_eq(1)./bs_vec_interp.^2;
% % neq_tau_pd = tau_array_neq(1)./bs_vec_interp.^2;
% 
% plot(ng_vec,tau_array_eq_ng,'Color',cmap2(3,:),'LineWidth',3)
% plot(ng_vec,tau_array_neq_ng,'Color',cmap2(2,:),'LineWidth',3)
% 
% % make scatters
% seq = scatter(ng_vec,tau_array_eq_ng,75,'MarkerFaceColor',cmap2(3,:),'MarkerEdgeColor','k');
% sneq = scatter(ng_vec,tau_array_neq_ng,75,'MarkerFaceColor',cmap2(2,:),'MarkerEdgeColor','k');
% 
%  ylim([30 1e4])
%  set(gca,'yscale','log')
% % legend([seq sneq],'equilibrium','non-equilibrium','Location','southwest')
% xlabel('number of general TF reactions');
% ylabel('decision time (burst cycles)');
% 
% grid on
% set(gca,'FontSize',14)
% set(gca,'Color',[228,221,209]/255) 
% 
% ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.XAxis(1).Color = 'k';
% 
% eq_neq_dtime_ng.InvertHardcopy = 'off';
% set(gcf,'color','w');
% 
% saveas(eq_neq_dtime_ng,[FigPath 'decision_time_vs_ng.png'])
% saveas(eq_neq_dtime_ng,[FigPath 'decision_time_vs_ng.pdf']) 
% 
% 
% %% Plot decision time for envelope
% 
% % initialize arrays to store rates
% tau_array_neq_g = NaN(1,length(phi_axis_ds_g));
% minError = 0.1;
% n_tries = 1e3;
% 
% % iterate
% for i = 1:length(phi_axis_g)
%        
%     % find best performer of 3 candidates
%     if any(~isnan(ir_gen_array(i,1:3)))
%         [~,mi] = nanmax(ir_gen_array(i,1:3));
% 
%         optimal_rates_neq = reshape(ir_gen_rates_array(i,mi,:),1,[]);
%         optimal_rates_neq = optimal_rates_neq(~isnan(optimal_rates_neq));
% 
%         % set function path
%         sweep_info_temp = master_struct_multi_g(mi).sweep_info;
%         functionPath = sweep_info_temp.functionPath;
%         slashes = strfind(functionPath,'\');
%         simName = functionPath(slashes(end-1)+1:slashes(end)-1);
%         rmpath(genpath('../utilities/metricFunctions/'));
%         addpath(genpath(['../utilities/metricFunctions/numeric/' simName]));
% 
%         sweep_info_temp.minError = minError;
%         sweep_info_temp.numerical_precision = 5;
% 
%         tau_neq = NaN;
%         iter = 1;
%         while isnan(tau_neq) && iter < n_tries
%             metric_vec_neq = calculateMetricsNumeric_v3(optimal_rates_neq, sweep_info_temp);
%             tau_neq = 1./metric_vec_neq(inv_dtime_index_num);
%             iter = iter + 1;
%         end
%         tau_array_neq_g(i) = tau_neq;    
%     end
% end
% 
% %% Make figure
% % fit smoothing spline
% [f_spline,~,~] = fit(log10(phi_axis_g'),tau_array_neq_g','smoothingspline','SmoothingParam',.85);
% f_spline_ds = feval(f_spline,log10(phi_axis_ds_g));
% 
% topline = repelem(1e3,length(f_spline_ds));
% 
% tau_g_fig = figure;
% hold on
% 
% % plot eq and noneq regimes
% fill([phi_axis_ds_g fliplr(phi_axis_ds_g)], [f_spline_ds' topline],cmap_pu(4,:),...
%                                               'FaceAlpha',0.5,'EdgeAlpha',0)
% plot(phi_axis_ds_g,f_spline_ds,'--k','LineWidth',4);  
% 
% set(gca,'xscale','log')
% ylim([0 750])
% xlim([5e-3 1e4])
% 
% xlabel('energy dissipation rate (k_BT per burst cycle)');
% ylabel('decision time (burst cycles)');
% grid on
% set(gca,'FontSize',14)
% set(gca,'Color',[228,221,209]/255) 
% 
% ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.XAxis(1).Color = 'k';
% 
% tau_g_fig.InvertHardcopy = 'off';
% set(gcf,'color','w');
% set(gca,'xscale','log')
% % set(gca,'yscale','log')
% 
% saveas(tau_g_fig,[FigPath 'dtime_vs_phi_g.png'])
% saveas(tau_g_fig,[FigPath 'dtime_vs_phi_g.pdf']) 
% 
% %%
% ir1 = nanmax(ir_gen_array(:,1))*0.98;
% phi_max = [];
% ir_max = [];
% for i = 1:size(ir_gen_array,2)
%     ir_max(i) = ir1*i;%max(ir_gen_array_sm(:,i));    
%     mi = find(ir_gen_array(:,i)>=ir1*i,1);
%     phi_max(i) = phi_axis_g(mi);
% end    