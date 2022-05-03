% script to conduct sweeps tracking IR vs pon for different CW values

clear 
close all
addpath(genpath('../utilities/'))

% set basic parameters
nStates = 6;
paramBounds = repmat([-10 ; 6],1,11); % constrain transition rate magnitude
[~,~,metric_names] = calculateMetricsSym_v2([]);

% specify function path
functionPath = ['../utilities/metricFunctions/n' num2str(nStates) '_OR/'];

% make sure we're linked to the appropriate function subfolder% make sure we're linked to the appropriate function subfolder
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));

% Set Dropbox directory
DropboxFolder = 'S:\Nick\Dropbox\Nonequilibrium\Nick\SweepOutput';
readPath = [DropboxFolder filesep 'sweeps06_r_vs_info_vs_cw' filesep];


% get index of useful metrics
ir_index = find(strcmp(metric_names,'DecisionRateNorm'));
sharp_index = find(strcmp(metric_names,'Sharpness'));
rate_index = find(strcmp(metric_names,'Production Rate'));
cw_index = find(strcmp(metric_names,'CW'));
rate_ent_index = find(strcmp(metric_names,'rateEntropy'));

% set sim options
sweep_options = {'n_seeds',5,'n_iters_max',50,'n_sim',5,'nStates',nStates};

% calculate sensitivity bound
alpha_factor = 100;

% specify plot options 
markerAlpha = 0.5; % marker transparency
markerSize = 75; % marker size
n_plot = 3e3;
bit_factor = log2(exp(1));

cw_vec = logspace(0,6,151);
n_keep = 1e2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read in parameter sweep results
rate_bins = 0.075:0.05:0.925;
rate_axis = rate_bins(1:end-1)+diff(rate_bins)/2;

% generate indexing vectors
rate_index_vec = repelem(1:length(rate_axis),length(cw_vec))';
cw_index_vec = repmat(1:length(cw_vec),1,length(rate_axis))';


rate_index_vec0 = repmat(1:length(rate_axis),1,length(cw_vec))';
cw_index_vec0 = repelem(1:length(cw_vec),length(rate_axis))';
rate_index_vec0 = repelem(rate_index_vec0,n_keep);
cw_index_vec0 = repelem(cw_index_vec0,n_keep);

rate_array_neq = NaN(length(cw_vec)*length(rate_axis),11);
rate_array_eq = NaN(length(cw_vec)*length(rate_axis),11);
rate_array_eq0 = NaN(length(cw_vec)*length(rate_axis)*n_keep,11);

ir_array_neq = NaN(length(cw_vec),length(rate_bins)-1);
ir_array_eq = NaN(length(cw_vec),length(rate_bins)-1);
ir_array_eq0 = NaN(length(cw_vec),length(rate_bins)-1,n_keep);

master_struct = struct;
try
  parpool(24)
catch  
  % do nothing
end  
  
parfor c = 1:length(cw_vec)
    
    suffix = ['cw' sprintf('%03d',c)];
    
    % load neq
    master_struct(c).sweep_info_neq = load([readPath 'sweep_info_neq_' suffix '.mat'], 'sweep_info_neq');
    master_struct(c).sweep_results_neq = load([readPath 'sweep_results_neq_' suffix '.mat'], 'sweep_results_neq');
    
    % load eq
    master_struct(c).sweep_info_eq = load([readPath 'sweep_info_eq_' suffix '.mat'], 'sweep_info_eq');
    master_struct(c).sweep_results_eq = load([readPath 'sweep_results_eq_' suffix '.mat'], 'sweep_results_eq');
    
    % load eq generic
    master_struct(c).sweep_info_eq0 = load([readPath 'sweep_info_eq0_' suffix '.mat'], 'sweep_info_eq0');
    master_struct(c).sweep_results_eq0 = load([readPath 'sweep_results_eq0_' suffix '.mat'], 'sweep_results_eq0');
        
end   
%%
first_i = 1;

for c = 1:length(cw_vec)    
    suffix = ['cw' sprintf('%03d',c)];
        
    sweep_results_neq = master_struct(c).sweep_results_neq.sweep_results_neq;    
    sweep_results_eq = master_struct(c).sweep_results_eq.sweep_results_eq;
    sweep_results_eq0 = master_struct(c).sweep_results_eq0.sweep_results_eq0;
    
    % find best gene circuits as a function of pon
    for r = 1:length(rate_bins)-1
        % find systems with correct range of production rates
        rate_vec_neq = sweep_results_neq.metric_array(:,rate_index);
        sharp_vec_neq = sweep_results_neq.metric_array(:,sharp_index);
        rate_filter_neq1 = rate_vec_neq>=rate_bins(r) & rate_vec_neq<rate_bins(r+1) & sharp_vec_neq>=0;
        ir_vec_neq = sweep_results_neq.metric_array(:,ir_index);
        max_ir_neq = nanmax(ir_vec_neq(rate_filter_neq1));
        rate_filter_neq = find(rate_filter_neq1 & ir_vec_neq <= 0.75*max_ir_neq);
        
        rate_vec_eq = sweep_results_eq.metric_array(:,rate_index);
        sharp_vec_eq = sweep_results_eq.metric_array(:,sharp_index);        
        rate_filter_eq = find(rate_vec_eq>=rate_bins(r) & rate_vec_eq<rate_bins(r+1) & sharp_vec_eq >= 0);
        
        rate_vec_eq0 = sweep_results_eq0.metric_array(:,rate_index);
        sharp_vec_eq0 = sweep_results_eq0.metric_array(:,sharp_index);       
        rate_filter_eq0 = find(rate_vec_eq0>=rate_bins(r) & rate_vec_eq0<rate_bins(r+1) & sharp_vec_eq0 >= 0);
    
        % find optimal IR systems 
        [ir_array_neq(c,r), max_i_neq] = nanmax(sweep_results_neq.metric_array(rate_filter_neq,ir_index));
        [ir_array_eq(c,r), max_i_eq] = nanmax(sweep_results_eq.metric_array(rate_filter_eq,ir_index));
        
        % randomly select generic eq systems
        eq0_ids = randsample(rate_filter_eq0,n_keep,true);
        ir_array_eq0(c,r,:) = reshape(sweep_results_eq0.metric_array(eq0_ids,ir_index),1,1,[]);
        
        % store system transition rates
        ind = (r-1)*length(cw_vec) + c;
        rate_array_neq(ind,:) = sweep_results_neq.rate_array(rate_filter_neq(max_i_neq),:);
        rate_array_eq(ind,:) = sweep_results_eq.rate_array(rate_filter_eq(max_i_eq),:);
        
        rate_array_eq0(first_i:first_i+n_keep-1,:) = sweep_results_eq0.rate_array(eq0_ids,:);
        first_i = first_i + n_keep;
    end    
end   
                                         
%% Generate mutated systems 
sweep_info_neq = master_struct(1).sweep_info_neq.sweep_info_neq;
sweep_info_eq0 = master_struct(1).sweep_info_eq0.sweep_info_eq0;

% add functions to path
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));

% set basic parameters
n_perturbations = 1e2;
m_factor = logspace(0,log10(alpha_factor),n_perturbations);
alpha_adjusted = alpha_factor ./ m_factor;

% generate arrays of optimal rates

% neq global
ir_rate_array_neq = repmat(rate_array_neq,1,1, length(m_factor));
ir_rate_array_neq(:,sweep_info_neq.unbindingFlags==1,:) = ir_rate_array_neq(:,sweep_info_neq.unbindingFlags==1,:) .* reshape(m_factor,1,1,[]);
ir_rate_array_neq(:,sweep_info_neq.b_index,:) = repmat(reshape(alpha_adjusted,1,1,[]),size(ir_rate_array_neq,1),1,1);

% now do this for the equilibrium "hord"
eq_rate_array = repmat(rate_array_eq0, 1,1,length(m_factor));
eq_rate_array(:,sweep_info_eq0.unbindingFlags==1,:) = eq_rate_array(:,sweep_info_eq0.unbindingFlags==1,:) .* reshape(m_factor,1,1,[]);
eq_rate_array(:,sweep_info_eq0.b_index,:) = repmat(reshape(alpha_adjusted,1,1,[]),size(eq_rate_array,1),1,1);


% calculate predicted observed sharpness and production for each kind of 
% network under different perturbation strengths

sharpness_vec_neq = NaN(size(ir_rate_array_neq,1),size(ir_rate_array_neq,3));
sharpness_vec_eq0 = NaN(size(eq_rate_array,1),size(eq_rate_array,3));

rate_vec_neq = NaN(size(ir_rate_array_neq,1),size(ir_rate_array_neq,3));
rate_vec_eq0 = NaN(size(eq_rate_array,1),size(eq_rate_array,3));

parfor p = 1:size(ir_rate_array_neq,3)
  
    % extract rates   
    neq_rates = ir_rate_array_neq(:,:,p);         
    eq_rates0 = eq_rate_array(:,:,p);   
    
    % convert to cell array and run calculations    
    paramCellNeq = mat2cell(neq_rates,size(neq_rates,1),ones(1,size(neq_rates,2)));    
    paramCellEq0 = mat2cell(eq_rates0,size(eq_rates0,1),ones(1,size(eq_rates0,2)));    
    
    % production rate     
    rate_vec_neq(:,p) = productionRateFunction(paramCellNeq{:});    
    rate_vec_eq0(:,p) = productionRateFunction(paramCellEq0{:});   
    
    % sharpness    
    sharpness_vec_neq(:,p) = sharpnessFunction(paramCellNeq{:});%./ (rate_vec_neq(:,p).*(1-rate_vec_neq(:,p)));        
    sharpness_vec_eq0(:,p) = sharpnessFunction(paramCellEq0{:});% ./ (rate_vec_eq_full(:,p).*(1-rate_vec_eq_full(:,p)));        
   
end

% Let's see what the predicted sharpness shift looks like for different
% designate index of perturbation strength to use for plots
close all

% normalize arrays
sharpness_array_norm_neq = sharpness_vec_neq ./ sharpness_vec_neq(:,1);
sharpness_array_norm_eq0 = sharpness_vec_eq0 ./ sharpness_vec_eq0(:,1);

rate_array_norm_neq = rate_vec_neq ./ rate_vec_neq(:,1);
rate_array_norm_eq0 = rate_vec_eq0 ./ rate_vec_eq0(:,1);

%%

% set plot filters
m_ind = 16;
r_ind = 9;

rate_filter0 = rate_index_vec0 == r_ind;%true(size(rate_index_vec0 >=r_ind-2 & rate_index_vec0 <= r_ind+2));
rate_filter = rate_index_vec == r_ind;



% close all
figure;
hold on
scatter(rate_array_norm_eq0(rate_filter0,m_ind),sharpness_array_norm_eq0(rate_filter0,m_ind).*m_factor(m_ind),50)
% scatter(rate_array_norm_eq0(rate_filter0,2*m_ind),sharpness_array_norm_eq0(rate_filter0,2*m_ind).*m_factor(2*m_ind),50)
% scatter(imgaussfilt(rate_array_norm_neq(rate_filter,m_ind),1), imgaussfilt(sharpness_array_norm_neq(rate_filter,m_ind),1).*m_factor(m_ind),50,'filled')
scatter(rate_array_norm_neq(rate_filter,m_ind), sharpness_array_norm_neq(rate_filter,m_ind).*m_factor(m_ind),50,'filled','s')
% scatter(rate_array_norm_neq(rate_filter,2*m_ind), sharpness_array_norm_neq(rate_filter,2*m_ind).*m_factor(2*m_ind),50,'filled','s')
%%
close all

sharp_shift_fig = figure;
cmap = brewermap([],'Set2');
cmap2 = flipud(brewermap([],'Spectral'));
colormap(cmap2)
hold on
scatter(cw_vec(cw_index_vec0(rate_filter0)),sharpness_array_norm_eq0(rate_filter0,m_ind).*m_factor(m_ind),'s','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k')
scatter(cw_vec(cw_index_vec(rate_filter)),imgaussfilt(sharpness_array_norm_neq(rate_filter,m_ind).*m_factor(m_ind),1),markerSize,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k')

set(gca,'xscale','log')
xlabel('relative non-cognate factor concentration (c / c)');
ylabel('normalized sharpness shift (s/s \times k/k)')
xlim([1e0 alpha_factor^3])
ylim([0 4])

set(gca,'xtick',[1 alpha_factor^1 alpha_factor^2 alpha_factor^3 alpha_factor^4]);%,'xticklabels',{'\alpha^0','\alpha^1','\alpha^2','\alpha^3','\alpha^4'})
xlim([1 1e5])
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
sharp_shift_fig.InvertHardcopy = 'off';
set(gcf,'color','w');


rate_shift_fig = figure;
cmap = brewermap([],'Set2');
cmap2 = flipud(brewermap([],'Spectral'));
colormap(cmap2)
hold on
scatter(cw_vec(cw_index_vec0(rate_filter0)),rate_array_norm_eq0(rate_filter0,m_ind),'s','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k')
scatter(cw_vec(cw_index_vec(rate_filter)),imgaussfilt(rate_array_norm_neq(rate_filter,m_ind),1),markerSize,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k')

set(gca,'xscale','log')
xlabel('relative non-cognate factor concentration (c / c)');
ylabel('normalized sharpness shift (s/s \times k/k)')
xlim([1e0 alpha_factor^3])
ylim([0 1.5])

set(gca,'xtick',[1 alpha_factor^1 alpha_factor^2 alpha_factor^3 alpha_factor^4]);%,'xticklabels',{'\alpha^0','\alpha^1','\alpha^2','\alpha^3','\alpha^4'})
xlim([1 1e5])
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
rate_shift_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
% saveas(sharp_shift_fig,[FigPath 'sharp_fold_vs_cw.png'])
% saveas(sharp_shift_fig,[FigPath 'sharp_fold_vs_cw.pdf'])
%%
sharp_shift_fig = figure;%('Position',[100 100 512 256]);
cmap = brewermap([],'Set2');
cmap2 = flipud(brewermap([],'Spectral'));
colormap(cmap2)
hold on
scatter(rate_array_norm_neq_full(:,m_ind),sharpness_array_norm_neq_full(:,m_ind).*m_factor(m_ind),[],ir_keep_values>=0.95,'filled','o','MarkerEdgeColor','k','MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0)
scatter(rate_array_norm_eq0(:,m_ind),sharpness_array_norm_eq0(:,m_ind).*m_factor(m_ind),'s','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k')
scatter(imgaussfilt(rate_array_norm_neq(:,m_ind),1),imgaussfilt(sharpness_array_norm_neq(:,m_ind).*m_factor(m_ind),1),markerSize,'MarkerFaceColor',cmap(8,:),'MarkerEdgeColor','k')


%%
close all
% rate shift
rate_shift_fig = figure;%('Position',[100 100 512 256]);
cmap = brewermap([],'Set2');
cmap2 = flipud(brewermap([],'Spectral'));
colormap(cmap2)
ir_ft = ir_keep_values>0.99;
hold on
scatter(cw_vec_neq(ir_keep_indices),rate_array_norm_neq_full(:,m_ind),[],ir_keep_values,'filled','o','MarkerEdgeColor','k','MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0)
% scatter(cw_vec_neq(ir_keep_indices(ir_ft)),rate_array_norm_neq_full(ir_ft,m_ind),[],'filled','o','MarkerEdgeColor','k','MarkerFaceAlpha',1,'MarkerEdgeAlpha',0)
% scatter(cw_vec_eq0(eq0_keep_indices),rate_array_norm_eq0(:,m_ind),'s','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k')
scatter(cw_vec_plot,imgaussfilt(rate_array_norm_neq(:,m_ind),1),markerSize,'MarkerFaceColor',cmap(8,:),'MarkerEdgeColor','k')

set(gca,'xscale','log')
xlabel('relative non-cognate factor concentration (c / c)');
ylabel('rate shift (r/r)')

set(gca,'xtick',[1 alpha_factor^1 alpha_factor^2 alpha_factor^3 alpha_factor^4]);%,'xticklabels',{'\alpha^0','\alpha^1','\alpha^2','\alpha^3','\alpha^4'})

xlim([1e0 alpha_factor^4])
set(gca,'FontSize',14)
ylim([.3 1.1])
xlim([1 1e5])
set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
% set(gca,'yscale','log')
rate_shift_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

% saveas(rate_shift_fig,[FigPath 'rate_fold_vs_cw.png'])
% saveas(rate_shift_fig,[FigPath 'rate_fold_vs_cw.pdf'])
%% Phase space plot
close all

m_ind_vec = fliplr([m_ind 36 51]);

% generate colormaps
c_cell = cell(1,3);
c_cell{1} = brewermap(length(cw_vec_plot),'Purples');
c_cell{2} = brewermap(length(cw_vec_plot),'Greens');
c_cell{3} = brewermap(length(cw_vec_plot),'Reds');

% classify eq dots
cw_vec_eq0_plot = cw_vec_eq0(eq0_keep_indices);
cw_eq0_indices = discretize(cw_vec_eq0_plot,cw_vec);

phase_space_fig = figure;
hold on
s = [];
for m = 1:length(m_ind_vec)
    m_ind_2 = m_ind_vec(m);
    cmap = c_cell{m};
    r_vec = imgaussfilt(rate_array_norm_neq(:,m_ind_2),1);
    s_vec = imgaussfilt(sharpness_array_norm_neq(:,m_ind_2),1);
    for c = 1:length(cw_vec)-1
        scatter(rate_array_norm_eq0(cw_eq0_indices==c,m_ind_2),sharpness_array_norm_eq0(cw_eq0_indices==c,m_ind_2).*m_factor(m_ind_2),[],...
              's','MarkerFaceColor',cmap(c,:),'MarkerEdgeColor','k','MarkerEdgeAlpha',.1)
            
        if c == 26
            s(m) = scatter(r_vec(c),s_vec(c).*m_factor(m_ind_2),15 + c^2/15,...
                      'MarkerFaceColor',cmap(c,:),'MarkerEdgeColor','k','MarkerEdgeAlpha',.2,'MarkerFaceAlpha',0.75);
        else
            scatter(r_vec(c),s_vec(c).*m_factor(m_ind_2),15 + c^2/15,...
                      'MarkerFaceColor',cmap(c,:),'MarkerEdgeColor','k','MarkerEdgeAlpha',.2,'MarkerFaceAlpha',0.75);
        end
        
        if mod(c,5)==0 && m == 1
            scatter(.1+c/500,.25,15 + c^2/15,...
                      'MarkerFaceColor',cmap(c,:),'MarkerEdgeColor','k','MarkerEdgeAlpha',.2,'MarkerFaceAlpha',0.5);
        end
    end
end
% set(gca,'xscale','log')
set(gca,'yscale','log')
ylim([10^-1.2 1e2])
xlabel('rate shift (r^*/r)')
ylabel('normalized sharpness shift (s^*/s \times k_-^*/k_-^0)')
set(gca,'FontSize',14)
legend(s,'k_-^*/k_-^0 = 2','k_-^*/k_-^0 = 5','k_-^*/k_-^0 = 10','Location','northwest','Fontsize',10)
set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
set(gca,'yscale','log')
phase_space_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(phase_space_fig,[FigPath 'phase_space_fig.png'])
saveas(phase_space_fig,[FigPath 'phase_space_fig.pdf'])

%%
close all
% sharpness fold shift vs perturbation strength
alpha_shift_fig = figure;
cmap = brewermap([],'Set2');
hold on
% scatter(m_factor / alpha_factor,fold_s_cw_neq_plot,markerSize,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k')
% scatter(m_factor / alpha_factor,fold_s_cw_eq_plot,markerSize,'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k')
plot(m_factor / alpha_factor,fold_s_cw_neq_plot,'Color',cmap(2,:),'LineWidth',3)
plot(m_factor / alpha_factor,fold_s_cw_eq_plot,'Color',cmap(3,:),'LineWidth',3)
set(gca,'yscale','log')
set(gca,'xscale','log')
xlabel('perturbation strength (\alpha^* / \alpha)');
ylabel('fold sharpness decrease (s^*/s)')
% xlim([1e0 alpha_factor^4])
ylim([10^-4 10^0])

% set(gca,'xtick',[1 alpha_factor^1 alpha_factor^2 alpha_factor^3 alpha_factor^4]);%,'xticklabels',{'\alpha^0','\alpha^1','\alpha^2','\alpha^3','\alpha^4'})

set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
alpha_shift_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
saveas(alpha_shift_fig,[FigPath 'sharp_fold_vs_alpha.png'])
saveas(alpha_shift_fig,[FigPath 'sharp_fold_vs_alpha.pdf'])
