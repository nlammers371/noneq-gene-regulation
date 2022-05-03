% Plot results for IR vs energy for higher order models
clear 
close all
addpath(genpath('../utilities/'))

% %%%%%%%%%%%%%%%%  set relative read and write paths %%%%%%%%%%%%%%%%%%%%
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\';
DataPath = [DropboxFolder  'SweepOutput\sweeps01_info_vs_energy' filesep ];
FigPath = [DropboxFolder '\presentatations' filesep 'UChicago' filesep];
mkdir(FigPath);

% %%%%%%%%%%%%%%%%  Set plot parameters and constants %%%%%%%%%%%%%%%%%%%%
n_plot = 3e3; % number of points to plot
% markerAlpha = 0.5; % marker transparency
% markerSize = 75; % marker size

close all
ds_factor = 3;
color_ind = 5;
sym_list = {'o','^','d','s','v'};
bs_vec = 1:length(master_struct_multi_bs);
t_cycle = 5;
% set plot parameters
markerSize = 50;
rng(231);
alphaFactor = 0.25;
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

% calculate upper IR bounds for each NBS
phi_bins = logspace(-3,log10(50),201);
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


% phi_max = 50; % max energy to plot in kBT
 
% make figure

% Define colormaps for use throughout
cmap_pu = brewermap(8,'Purples');
cmap_rd = brewermap(8,'Reds');
cmap_bu = brewermap(8,'Blues');
cmap_gre = brewermap(8,'Greens');
cmap_gra = brewermap(8,'Greys');

close all

cmap = [cmap_pu(color_ind,:); cmap_gre(color_ind,:); cmap_rd(color_ind,:); ...
          cmap_bu(color_ind,:) ; cmap_gra(color_ind,:)];

         
% generate down-sampled version of x axis for plotting
phi_bs_ds = downsample(phi_axis,ds_factor);
ir_bound_array_fit = NaN(length(phi_bs_ds),size(ir_bound_array,2));

% fit smoothing splines to upper bounds for each genecircuit type

for i = length(master_struct_multi_bs):-1:1  
    % fit smoothing spline
    ir_bound_temp = ir_bound_array(:,i);
    ir_bound_temp = ir_bound_temp(~isnan(ir_bound_temp));
    phi_temp = phi_axis(~isnan(ir_bound_array(:,i)));
    phi_max = max(phi_temp);
    phi_bs_ds_temp = phi_bs_ds(phi_bs_ds<=phi_max);
    
    [f_spline,~,~] = fit(log10(phi_temp'),ir_bound_temp,'smoothingspline','SmoothingParam',0.99);
    f_spline_ds = feval(f_spline,log10(phi_bs_ds_temp));%downsample(f_spline,ds_factor);    

    % plot smoothed trends
    ir_bound_array_fit(phi_bs_ds<=phi_max,i) = f_spline_ds; 
end    

% remove phi values for which one more more IR values are missing
nan_flags = any(isnan(ir_bound_array_fit),2);
ir_bound_array_fit = ir_bound_array_fit(~nan_flags,:);
phi_bs_ds = phi_bs_ds(~nan_flags);

% %%%%%%%%%%%% Look at impact of multiple General TFs %%%%%%%%%%%%%%%%%%%%

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

%%
% Make plots for presentation
close all

[~,eq_ind_g] = min(abs(phi_axis_ds_g-0.005));
[~,eq_ind_bs] = min(abs(phi_bs_ds-0.005));



% "equilibrium" vs. "non-equilibrium" strategies
eq_neq_fig = figure;
hold on

% ir1 = ir_gen_array_sm(eq_ind_g,1)/t_cycle;
% phi1 = phi_axis_ds_g(eq_ind_g)/t_cycle;
% scatter(phi1,ir1,75,'MarkerFaceColor',brighten(cmap_pu(3,:),-0.25),'MarkerEdgeColor','k'); 

set(gca,'xscale','log')
ylim([0 0.07]/t_cycle)
xlim([4e-3 1e4]/t_cycle)

xlabel('energy dissipation rate (k_BT per minute)');
ylabel('information rate (bits per minute)')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
set(gcf,'color','w');
eq_neq_fig.InvertHardcopy = 'off';

s = [];
% layer on neq LC models
for i = 1:size(ir_gen_array_sm,2)
  
    ir_bound = ir_gen_array_sm(:,i)/t_cycle;
    nan_filter = ~isnan(ir_bound)&(1:length(phi_axis_ds_g)>=eq_ind_g)';
    phi_plot = phi_axis_ds_g(nan_filter)/t_cycle;
    ir_plot = ir_bound(nan_filter);
    
    fill([phi_plot fliplr(phi_plot)], [ir_plot' repelem(1e-6,length(ir_plot))],cmap_pu(2+i,:),...
                                        'FaceAlpha',1,'EdgeAlpha',0);
                                      
    s(end+1) = scatter(phi_plot(1:2:end),ir_plot(1:2:end),75,'MarkerFaceColor',brighten(cmap_pu(2+i,:),-0.25),'MarkerEdgeColor','k');
    
    h = get(gca,'Children');
    set(gca,'Children',[h(3:end)' h(1) h(2)])

    if i == 4
        legend(s,'N_g=2','N_g=3','N_g=4','N_g=5','Location','northwest','Color','w');
    end
    
    save_name = ['ir_eq_neq_' sprintf('%02d',i)];
    saveas(eq_neq_fig,[FigPath save_name '.png'])
    saveas(eq_neq_fig,[FigPath save_name '.pdf']) 
        
end  

% add IR bound
ir_eq = nanmax(ir_gen_array(phi_axis_g<5e-3,1));
ir_bound_pd = ir_eq + ir_eq*sqrt(phi_axis_g);
plot(phi_axis_g/t_cycle,ir_bound_pd/t_cycle,':k','LineWidth',4);
legend(s,'N_g=2','N_g=3','N_g=4','N_g=5','Location','northwest','Color','w');
save_name = ['ir_neq_LC' sprintf('%02d',i+1)];
saveas(eq_neq_fig,[FigPath save_name '.png'])
saveas(eq_neq_fig,[FigPath save_name '.pdf']) 

%%
% "equilibrium" vs. "non-equilibrium" strategies

eq_neq_fig = figure;
hold on

set(gca,'xscale','log')
set(gca,'yscale','log')
ylim([1e-3 0.5]/t_cycle)
xlim([4e-3 3e1]/t_cycle)

xlabel('energy dissipation rate (k_BT per minute)');
ylabel('information rate (bits per minute)')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
set(gcf,'color','w');
eq_neq_fig.InvertHardcopy = 'off';

s = [];
% layer on eq bs models
for i = 1:length(bs_vec)
%     scatter(phi_bs_ds(eq_ind_bs)/t_cycle,ir_bound_array_fit(:,i)/t_cycle,75,'MarkerFaceColor',brighten(cmap(i,:),-0.25),'MarkerEdgeColor','k');

    % downsample vectors
    if i > 1
        fill([phi_bs_ds fliplr(phi_bs_ds)]/t_cycle,[ir_bound_array_fit(:,i)' fliplr(ir_bound_array_fit(:,i-1)')]/t_cycle,...
                                        cmap(i,:),'FaceAlpha',alphaFactor*2,'EdgeAlpha',0);
    else
        fill([phi_bs_ds fliplr(phi_bs_ds)]/t_cycle,[ir_bound_array_fit(:,i)' repelem(1e-6,length(ir_bound_array_fit(:,i)))]/t_cycle,...
                                        cmap(i,:),'FaceAlpha',alphaFactor*2,'EdgeAlpha',0);
    end


    % overlay line plots
    s(end+1) = scatter(phi_bs_ds(1:2:end)/t_cycle,ir_bound_array_fit(1:2:end,i)/t_cycle,markerSize','MarkerFaceColor',brighten(cmap(i,:),-0.25),'MarkerEdgeColor','k');
  
    h = get(gca,'Children');
    set(gca,'Children',[h(3:end)' h(1) h(2)])
    
    if i == length(bs_vec)
        legend(s,'N_a=1','N_a=2','N_a=3','N_a=4','N_a=5','Location','south','Orientation','horizontal','Color','w');
    end
    save_name = ['ir_eq_bs_' sprintf('%02d',i)];
    saveas(eq_neq_fig,[FigPath save_name '.png'])
    saveas(eq_neq_fig,[FigPath save_name '.pdf']) 
end    

%% plot corresponding decision times
close all
phi_axis = phi_bs_ds/t_cycle;

ir_bs = ir_bound_array_fit/t_cycle;
minError = 0.32;
e_factor = log((1-minError)./minError).*(1-2*minError);

dt_bs = e_factor./ir_bs;

topline = repelem(1e8,length(phi_axis));


dt_fig = figure;
hold on

ylim([1 600])
xlim([4e-3 3e1]/t_cycle)
set(gca,'yscale','log')
xlabel('energy dissipation rate (k_BT per minute)');
ylabel('decision time (minutes)')
grid on
set(gca,'FontSize',14)
set(gca,'xscale','log')

set(gca,'Color',[228,221,209]/255) 

set(gcf,'color','w');
dt_fig.InvertHardcopy = 'off';

for i = 1:length(bs_vec)
    
    if i > 1
        fill([phi_axis fliplr(phi_axis)],[dt_eq(:,i)' fliplr(dt_eq(:,i-1)')],...
                                        cmap(i,:),'FaceAlpha',alphaFactor*2,'EdgeAlpha',0);
    else
       fill([phi_axis fliplr(phi_axis)],[dt_eq(:,i)' topline],...
                                        cmap(i,:),'FaceAlpha',alphaFactor*2,'EdgeAlpha',0);
    end
    
    plot(phi_axis,dt_eq(:,i),'-','Color',brighten(cmap(i,:),-.2),'LineWidth',3); 
    
    h = get(gca,'Children');
    set(gca,'Children',[h(3:end)' h(1) h(2)])
    
    saveas(dt_fig,[FigPath 'dt_bs_phi' num2str(i) '.png'])
    saveas(dt_fig,[FigPath 'dt_bs_phi' num2str(i) '.pdf']) 
end


%% Now for LC
close all

ir_lc = nanmax(ir_gen_array_sm,[],2)/t_cycle;
phi_plot = phi_axis_ds_g(~isnan(ir_lc))/t_cycle;

dt_lc_neq = e_factor./ir_lc;
topline = repelem(1e8,length(phi_plot));

% make plot
dt_fig = figure;
cmap2 = brewermap([],'Set2');
hold on

ylim([1 600])
xlim([4e-3 1e4]/t_cycle)

set(gca,'yscale','log')
xlabel('energy dissipation rate (k_BT per minute)');
ylabel('decision time (minutes)')
grid on
set(gca,'FontSize',14)
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'Color',[228,221,209]/255) 

% set(gca,'xscale','log')
set(gcf,'color','w');
dt_fig.InvertHardcopy = 'off';

for i = 1%:size(dt_lc_neq,2)
    
    if i > 1
        fill([phi_plot fliplr(phi_plot)],[dt_lc_neq(:,i)' fliplr(dt_lc_neq(:,i-1)')],...
                                        cmap_pu(5,:),'FaceAlpha',alphaFactor*2,'EdgeAlpha',0);
    else
       fill([phi_plot fliplr(phi_plot)],[dt_lc_neq(:,i)' topline],...
                                        cmap_pu(5,:),'FaceAlpha',alphaFactor*3,'EdgeAlpha',0);
    end
    
    plot(phi_plot,dt_lc_neq(:,i),'-','Color',brighten(cmap_pu(5,:),-.2),'LineWidth',3); 
    
    h = get(gca,'Children');
    set(gca,'Children',[h(3:end)' h(1) h(2)])
    
    saveas(dt_fig,[FigPath 'dt_lc_neq' num2str(i) '.png'])
    saveas(dt_fig,[FigPath 'dt_lc_neq' num2str(i) '.pdf']) 
end

% plot IR PHI bound
dt_bound_pd = e_factor./ir_bound_pd *t_cycle;
plot(phi_axis_g/t_cycle,dt_bound_pd,':k','LineWidth',4);

saveas(dt_fig,[FigPath 'dt_lc_neq' num2str(i+1) '.png'])
saveas(dt_fig,[FigPath 'dt_lc_neq' num2str(i+1) '.pdf']) 
 