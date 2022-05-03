% Plot results for IR vs energy for higher order models
clear 
close all
addpath(genpath('../utilities/'))

% %%%%%%%%%%%%%%%%  set relative read and write paths %%%%%%%%%%%%%%%%%%%%
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\';
DataPath = [DropboxFolder  'SweepOutput\sweeps01_info_vs_energy' filesep ];
FigPath = [DropboxFolder '\manuscript\info_vs_energy' filesep];
mkdir(FigPath);

% %%%%%%%%%%%%%%%%  Set plot parameters and constants %%%%%%%%%%%%%%%%%%%%
n_plot = 3e3; % number of points to plot
markerAlpha = 0.5; % marker transparency
markerSize = 75; % marker size

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%% Examine impact of adding binding sites %%%%%%%%%%%%%%%%%
% get metric names for numeric sweeps
[~,~,metric_names_num] = calculateMetricsNumeric_v3([]);

ir_index_num = find(strcmp(metric_names_num,'IR'));
cycle_time_index_num = find(strcmp(metric_names_num,'TauCycle'));
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
phi_bins = logspace(-2.3,log10(50),201);
phi_axis = phi_bins(2:end);
ir_bound_array = NaN(length(phi_axis),length(master_struct_multi_bs));
ir_rates_array = NaN(length(phi_axis),length(master_struct_multi_bs),size(master_struct_multi_bs(end).sweep_results.rate_array,2));

for i = length(master_struct_multi_bs):-1:1
    metric_array = master_struct_multi_bs(i).sweep_results.metric_array;
    rate_array = master_struct_multi_bs(i).sweep_results.rate_array;
    phi_vec = metric_array(:,phi_index_num);
    ir_vec = metric_array(:,ir_index_num);
    
    for p = 1:length(phi_bins)-1
      
        phi_filter = phi_vec>phi_bins(p) & phi_vec<=phi_bins(p+1);
        ir_filtered = ir_vec(phi_filter);
        rates_filtered = rate_array(phi_filter,:);
        
        % deal with outliers due to numerical precision issues (NL: need to
        % resolve this)
        if phi_bins(p+1) <= 0.1
           ir_filtered(ir_filtered > 0.067) = NaN;
        end
        if any(ir_filtered)
            [ir_bound_array(p,i), mi] = nanmax(ir_filtered*log2(exp(1))); % convert to bits
            ir_rates = rates_filtered(mi,:);
            ir_rates_array(p,i,1:length(ir_rates)) = ir_rates;
        end
    end
end
%%
close all
ds_factor = 3;
color_ind = 5;
sym_list = {'o','^','d','s','v'};

% set plot parameters
markerSize = 50;
rng(231);
alphaFactor = 0.25;
phi_max = 33;
 
% make figure
info_bs_fig = figure;

% Define colormaps for use throughout
cmap_pu = brewermap(8,'Purples');
cmap_rd = brewermap(8,'Reds');
cmap_bu = brewermap(8,'Blues');
cmap_gre = brewermap(8,'Greens');
cmap_gra = brewermap(8,'Greys');

cmap = [cmap_pu(color_ind,:); cmap_gre(color_ind,:); cmap_rd(color_ind,:); ...
          cmap_bu(color_ind,:) ; cmap_gra(color_ind,:)];

        

 
% generate down-sampled version of x axis for plotting
phi_bs_ds = downsample(phi_axis,ds_factor);
ir_bound_array_sm = NaN(length(phi_bs_ds),size(ir_bound_array,2));
s = [];
hold on
for i = length(master_struct_multi_bs):-1:1    
  
    % extract raw vectors
    phi_vec_ds = master_struct_multi_bs(i).sweep_results.metric_array(:,phi_index_num);
    ir_vec_ds = master_struct_multi_bs(i).sweep_results.metric_array(:,ir_index_num)*log2(exp(1));
    
    % check for outliers
    out_flags = (phi_vec_ds<=0.1 & ir_vec_ds>=0.1)|phi_vec_ds<=5e-3|phi_vec_ds>=31|ir_vec_ds<0;
    use_indices = find(~isnan(ir_vec_ds)&~out_flags);
    
    % find bounary points    
    ds_factor = 3e3;
    grid_res = 3e2;
    
    [boundaryPoints, metric_array_filtered, random_grid] = findBoundaryPoints(...
                          [log10(phi_vec_ds) ir_vec_ds],use_indices,10,grid_res,ds_factor,1);
                        
    phi_vec_ds = phi_vec_ds(use_indices);
    ir_vec_ds = ir_vec_ds(use_indices);
    
    % fit smoothing spline
    ir_bound_temp = ir_bound_array(:,i);
    ir_bound_temp = ir_bound_temp(phi_axis<phi_max);
    phi_temp = phi_axis(phi_axis<phi_max);
   
    phi_bs_ds_temp = phi_bs_ds(phi_bs_ds<=phi_max);
    
    [f_spline,~,~] = fit(log10(phi_temp'),ir_bound_temp,'smoothingspline','SmoothingParam',0.95);
    f_spline_ds = feval(f_spline,log10(phi_bs_ds_temp));%downsample(f_spline,ds_factor);           
    
    plot_indices = [boundaryPoints ; random_grid];    
    plot_indices = plot_indices(~isnan(plot_indices));
        
    
    % downsample vectors
    scatter(phi_vec_ds(plot_indices),ir_vec_ds(plot_indices),...
                                        markerSize,sym_list{i},'MarkerFaceColor',cmap(i,:),...
                                        'MarkerEdgeColor','k','MarkerFaceAlpha',alphaFactor,...
                                        'MarkerEdgeAlpha',0.1);
                                      
    % plot smoothed trends
    ir_bound_array_sm(phi_bs_ds<=phi_max,i) = f_spline_ds;    
end

for i = 1:length(bs_vec)
    s(i) = plot(phi_bs_ds,ir_bound_array_sm(:,i),'Color',brighten(cmap(i,:),-0.25),'LineWidth',3);
end    

set(gca,'xscale','log')
ylim([1e-3 0.5])
xlim([1e-2 3e1])

% h = colorbar;
% ylabel(h,'number of binding sites')
legend(s,'1 bs','2 bs','3 bs','4 bs','5 bs','Location','southeast','Orientation','horizontal');
xlabel('energy dissipation rate (k_BT per burst cycle)');
ylabel('information rate (bits per burst cycle)')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

info_bs_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
set(gca,'xscale','log')
set(gca,'yscale','log')

saveas(info_bs_fig,[FigPath 'info_vs_phi_bs.png'])
saveas(info_bs_fig,[FigPath 'info_vs_phi_bs.pdf']) 

%% Plot eq and neq bounds bs ns
% note that I'm using low Phi as a proxy for eq performance
[~,eq_ind] = min(abs(phi_bs_ds-1e-2));
[~,neq_ind] = min(abs(phi_bs_ds-3e1));

eq_neq_bound_bs = figure;
cmap2 = brewermap(8,'Set2');
bs_vec = 1:length(master_struct_multi_bs);
hold on

ir_eq = ir_bound_array_sm(eq_ind,:);
ir_neq = ir_bound_array_sm(neq_ind,:);

% make plots
plot(bs_vec,ir_eq,'Color','k','LineWidth',1)
plot(bs_vec,ir_neq,'Color','k','LineWidth',1)

% make scatters
seq = scatter(bs_vec,ir_eq,75,'MarkerFaceColor',cmap2(3,:),'MarkerEdgeColor','k');
sneq = scatter(bs_vec,ir_neq,75,'MarkerFaceColor',cmap2(2,:),'MarkerEdgeColor','k');

% ylim([1e-3 0.5])
legend([seq sneq],'equilibrium','non-equilibrium','Location','northwest')
xlabel('number of binding sites');
ylabel('information rate (bits per burst cycle)')

grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

info_bs_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(info_bs_fig,[FigPath 'info_vs_bs.png'])
saveas(info_bs_fig,[FigPath 'info_vs_bs.pdf']) 

%% Do the same for decision times
% initialize arrays to store rates
% rate_array = NaN(n_lines,size(rate_array_flux,2));
[~,eq_ind_full] = min(abs(phi_axis-1e-2));
[~,neq_ind_full] = min(abs(phi_axis-3e1));

tau_array_neq = NaN(1,length(bs_vec));
tau_array_eq = NaN(1,length(bs_vec));
minError = 0.1;
n_tries = 1e3;
% iterate
for i = 1:length(bs_vec)
       
  optimal_rates_neq = reshape(ir_rates_array(neq_ind_full,i,:),1,[]);
  optimal_rates_neq = optimal_rates_neq(~isnan(optimal_rates_neq));
  optimal_rates_eq = reshape(ir_rates_array(eq_ind_full,i,:),1,[]);
  optimal_rates_eq = optimal_rates_eq(~isnan(optimal_rates_eq));

  % set function path
  sweep_info_temp = master_struct_multi_bs(i).sweep_info;
  functionPath = sweep_info_temp.functionPath;
  slashes = strfind(functionPath,'\');
  simName = functionPath(slashes(end-1)+1:slashes(end)-1);
  rmpath(genpath('../utilities/metricFunctions/'));
  addpath(genpath(['../utilities/metricFunctions/numeric/' simName]));

  sweep_info_temp.minError = minError;
  sweep_info_temp.numerical_precision = 5;
  
  tau_neq = NaN;
  iter = 1;
  while isnan(tau_neq) && iter < n_tries
      metric_vec_neq = calculateMetricsNumeric_v3(optimal_rates_neq, sweep_info_temp);
      tau_neq = 1./metric_vec_neq(inv_dtime_index_num);
      iter = iter + 1;
  end
  tau_array_neq(i) = tau_neq;
  
  tau_eq = NaN;
  iter = 1;
  while isnan(tau_eq) && iter < n_tries
      metric_vec_eq = calculateMetricsNumeric_v3(optimal_rates_eq, sweep_info_temp);
      tau_eq = 1./metric_vec_eq(inv_dtime_index_num);
      iter = iter + 1;
  end
  tau_array_eq(i) = tau_eq;
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make decision time vs binding site figure

topline = repelem(1e4,length(bs_vec));

eq_neq_dtime_bs = figure;
cmap2 = brewermap(8,'Set2');

hold on

% plot eq and noneq regimes
fill([bs_vec fliplr(bs_vec)], [tau_array_eq topline],cmap2(3,:),'FaceAlpha',0.5,'EdgeAlpha',0)
fill([bs_vec fliplr(bs_vec)], [tau_array_neq fliplr(tau_array_eq)],cmap2(2,:),'FaceAlpha',0.5,'EdgeAlpha',0)

% make plots
plot(bs_vec,tau_array_eq,'Color',cmap2(3,:),'LineWidth',3)
plot(bs_vec,tau_array_neq,'Color',cmap2(2,:),'LineWidth',3)

% make scatters
seq = scatter(bs_vec,tau_array_eq,75,'MarkerFaceColor',cmap2(3,:),'MarkerEdgeColor','k');
sneq = scatter(bs_vec,tau_array_neq,75,'MarkerFaceColor',cmap2(2,:),'MarkerEdgeColor','k');

 ylim([1e1 1e3])
 set(gca,'yscale','log')
% legend([seq sneq],'equilibrium','non-equilibrium','Location','southwest')
xlabel('number of binding sites');
ylabel('decision time (burst cycles)');

grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

eq_neq_dtime_bs.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(eq_neq_dtime_bs,[FigPath 'decision_time_vs_bs.png'])
saveas(eq_neq_dtime_bs,[FigPath 'decision_time_vs_bs.pdf']) 

%% Plot info gain vs Phi
info_gain_bs_fig = figure;

colormap(cmap(4:end,:));
hold on
for i = 1:length(master_struct_multi_bs)    
    % plot smoothed trends
    ir_sm = imgaussfilt(ir_bound_array(:,i),4);
%     plot(phi_axis,ir_sm/ir_sm(1),'Color',cmap(3+i,:),'LineWidth',3);
    scatter(phi_bs_ds,ir_bound_array_sm(:,i)/ir_bound_array_sm(eq_ind,i),50,sym_list{1},'MarkerFaceColor',cmap(i,:),...
              'MarkerEdgeColor','k');
end

set(gca,'xscale','log')
ylim([0 5])
xlim([1e-2 3e1])

% h = colorbar;
% ylabel(h,'number of binding sites')
legend('1 site','2 sites','3 sites','4 sites','5 sites','Location','northwest')
xlabel('energy dissipation rate (k_BT per burst cycle)');
ylabel('information rate (bits per burst cycle)')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

info_gain_bs_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
set(gca,'xscale','log')

saveas(info_gain_bs_fig,[FigPath 'app_info_gain_vs_bs.png'])
saveas(info_gain_bs_fig,[FigPath 'app_info_gain_vs_bs.pdf']) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%% Look at impact of multiple General TFs %%%%%%%%%%%%%%%%%%%%

% get list of sweep results files with only 1 genera TF reaction
multi_g_sweep_files = dir([DataPath 'sweep_results_s01_ns00_g0*']);
multi_g_sweep_files = [multi_g_sweep_files ; dir([DataPath 'sweep_results_s02_ns00_g0*'])];
multi_g_info_files = dir([DataPath 'sweep_info_s01_ns00_g0*']);
multi_g_info_files = [multi_g_info_files ; dir([DataPath 'sweep_info_s02_ns00_g0*'])];
  
% load
master_struct_multi_g = struct;
for f = 1:length(multi_g_sweep_files)
  
    load([DataPath multi_g_sweep_files(f).name])
    load([DataPath multi_g_info_files(f).name])
    
    master_struct_multi_g(f).sweep_results = sim_results;
    master_struct_multi_g(f).sweep_info = sim_info;
end    

%% calculate upper IR bounds

% calculate upper IR bounds for each NBS
phi_bins_g = logspace(-2.3,log10(1e4),51);
phi_axis_g = phi_bins_g(2:end);
ir_gen_array = NaN(length(phi_axis_g),length(master_struct_multi_g));
ir_gen_rates_array = NaN(length(phi_axis_g),length(master_struct_multi_g),size(master_struct_multi_g(end).sweep_results.rate_array,2));

for i = length(master_struct_multi_g):-1:1
    metric_array = master_struct_multi_g(i).sweep_results.metric_array;
    rate_array = master_struct_multi_g(i).sweep_results.rate_array;
    phi_vec = metric_array(:,phi_index_num);
    ir_vec = metric_array(:,ir_index_num);
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

%% Make figure 

% designate "cap" Phi values for different numbers of general reactions
phi_caps = [40 7e3 7e3];
ds_factor_3 = 1;

% downsample phi axis
phi_axis_ds_g = downsample(phi_axis_g,ds_factor_g);
ir_gen_array_sm = NaN(length(phi_axis_ds_g),size(ir_gen_array,2));

info_g_fig = figure;
s = [];
colormap(cmap(4:end,:));
hold on
for i = 3:-1:1    
   
    c_factor = color_ind;% + inc_factor;
    
    % fit spline to generate smoothed trend
    ir_bound_raw = ir_gen_array(:,i);        
    phi_filter = phi_axis_g<=phi_caps(i);
    phi_axis_ds_temp = phi_axis_ds_g(phi_axis_ds_g<=phi_caps(i));
    
    [f_spline,~,~] = fit(log10(phi_axis_g(phi_filter)'),ir_bound_raw(phi_filter),'smoothingspline','SmoothingParam',0.999);
    f_spline_ds = feval(f_spline,log10(phi_axis_ds_temp));
    
    phi_vec_ds = master_struct_multi_g(i).sweep_results.metric_array(:,phi_index_num);
    ir_vec_ds = master_struct_multi_g(i).sweep_results.metric_array(:,ir_index_num)*log2(exp(1));
    use_indices = find(~isnan(ir_vec_ds)&phi_vec_ds<=phi_caps(i));
    
    [boundaryPoints, metric_array_filtered, random_grid] = findBoundaryPoints(...
                          [log10(phi_vec_ds) ir_vec_ds],use_indices,10,3e2,3e3,1);
                        
                        
    plot_indices = [boundaryPoints ; random_grid];    
    plot_indices = plot_indices(~isnan(plot_indices));
    
    phi_vec_ds = phi_vec_ds(use_indices(plot_indices));
    ir_vec_ds = ir_vec_ds(use_indices(plot_indices));
    
    % plot smoothed trends
    ir_gen_array_sm(phi_axis_ds_g<=phi_caps(i),i) = f_spline_ds;
    
    % downsample vectors
    scatter(phi_vec_ds,ir_vec_ds,markerSize,'MarkerFaceColor',cmap_pu(1+i*2,:),...
                                        'MarkerEdgeColor','k','MarkerFaceAlpha',alphaFactor,...
                                        'MarkerEdgeAlpha',0.1);
                                      
    s(i) = plot(phi_axis_ds_temp,f_spline_ds,'Color',brighten(cmap_pu(1+i*2,:),-0.5),'LineWidth',3);    

%     plot(phi_axis_g,ir_gen_array_sm(:,i),'Color',cmap(ind,:),'LineWidth',3);
end

% fit spline to full set to obtain master curve
ir_gen_array_max = nanmax(ir_gen_array(:,1:3),[],2);
[f_spline,~,~] = fit(log10(phi_axis_g'),ir_gen_array_max,'smoothingspline','SmoothingParam',.99);
f_spline_ds = feval(f_spline,log10(phi_axis_ds_g));
plot(phi_axis_ds_g,f_spline_ds,'--k','LineWidth',4);  

set(gca,'xscale','log')
ylim([0 0.05])
xlim([1e-2 5e3])

% h = colorbar;
% ylabel(h,'number of binding sites')
legend(s,'1 general TF','2 general TFs','3 general TFs','Location','northwest');
xlabel('energy dissipation rate (k_BT per burst cycle)');
ylabel('information rate (bits per burst cycle)')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

info_g_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
set(gca,'xscale','log')
% set(gca,'yscale','log')

saveas(info_g_fig,[FigPath 'info_vs_phi_g.png'])
saveas(info_g_fig,[FigPath 'info_vs_phi_g.pdf']) 

%% Plot decision time for envelope

% initialize arrays to store rates
tau_array_neq_g = NaN(1,length(phi_axis_ds_g));
minError = 0.1;
n_tries = 1e3;

% iterate
for i = 1:length(phi_axis_g)
       
    % find best performer of 3 candidates
    if any(~isnan(ir_gen_array(i,1:3)))
        [~,mi] = nanmax(ir_gen_array(i,1:3));

        optimal_rates_neq = reshape(ir_gen_rates_array(i,mi,:),1,[]);
        optimal_rates_neq = optimal_rates_neq(~isnan(optimal_rates_neq));

        % set function path
        sweep_info_temp = master_struct_multi_g(mi).sweep_info;
        functionPath = sweep_info_temp.functionPath;
        slashes = strfind(functionPath,'\');
        simName = functionPath(slashes(end-1)+1:slashes(end)-1);
        rmpath(genpath('../utilities/metricFunctions/'));
        addpath(genpath(['../utilities/metricFunctions/numeric/' simName]));

        sweep_info_temp.minError = minError;
        sweep_info_temp.numerical_precision = 5;

        tau_neq = NaN;
        iter = 1;
        while isnan(tau_neq) && iter < n_tries
            metric_vec_neq = calculateMetricsNumeric_v3(optimal_rates_neq, sweep_info_temp);
            tau_neq = 1./metric_vec_neq(inv_dtime_index_num);
            iter = iter + 1;
        end
        tau_array_neq_g(i) = tau_neq;    
    end
end

%% Make figure
% fit smoothing spline
[f_spline,~,~] = fit(log10(phi_axis_g'),tau_array_neq_g','smoothingspline','SmoothingParam',.99);
f_spline_ds = feval(f_spline,log10(phi_axis_ds_g));

topline = repelem(1e3,length(f_spline_ds));

tau_g_fig = figure;
hold on

% plot eq and noneq regimes
fill([phi_axis_ds_g fliplr(phi_axis_ds_g)], [f_spline_ds' topline],cmap_pu(4,:),...
                                              'FaceAlpha',0.5,'EdgeAlpha',0)
plot(phi_axis_ds_g,f_spline_ds,'--k','LineWidth',4);  

set(gca,'xscale','log')
ylim([0 700])
xlim([1e-2 5e3])

xlabel('energy dissipation rate (k_BT per burst cycle)');
ylabel('decision time (burst cycles)');
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

tau_g_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
set(gca,'xscale','log')
% set(gca,'yscale','log')

saveas(tau_g_fig,[FigPath 'dtime_vs_phi_g.png'])
saveas(tau_g_fig,[FigPath 'dtime_vs_phi_g.pdf']) 
