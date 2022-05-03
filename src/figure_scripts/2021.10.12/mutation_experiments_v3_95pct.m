% script to track information rate as a function of cw

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

% define save path

DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\manuscript\';
% DropboxFolder = 'S:\Nick\Dropbox\Nonequilibrium\Nick\manuscript\';

FigPath = [DropboxFolder 'observed_sharpness' filesep];
mkdir(FigPath);         

% get index of useful metrics
spec_index = find(strcmp(metric_names,'Specificity'));
sharp_right_norm_index = find(strcmp(metric_names,'SharpnessRightNorm'));
ir_index = find(strcmp(metric_names,'DecisionRateNorm'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
sharpness_norm_index = find(strcmp(metric_names,'SharpnessNormed'));
precision_index = find(strcmp(metric_names,'Precision'));
rate_index = find(strcmp(metric_names,'Production Rate'));
precision_right_index = find(strcmp(metric_names,'PrecisionRight'));
cw_index = find(strcmp(metric_names,'CW'));
rate_ent_index = find(strcmp(metric_names,'rateEntropy'));

% set sim options
sweep_options = {'n_seeds',5,'n_iters_max',50,'n_sim',10,'nStates',nStates};

% calculate sensitivity bound
alpha_factor = 100;

% specify plot options 
markerAlpha = 0.5; % marker transparency
markerSize = 75; % marker size
n_plot = 3e3;
bit_factor = log2(exp(1));

% Load data
% load([FigPath 'mutation_experiment_data'],'master_struct')

%% Simulate perturbations to binding site of differing levels of severity
%%%%%%%%%%%%%%%%%%%%
% conduct sweeps
tic
% global neq
[sim_info_neq, sim_struct_neq] = param_sweep_multi_v3([cw_index ir_index],functionPath,sweep_options{:},...
                                              'half_max_flag',true,...'cw',1,...
                                              'equilibrium_flag',false,'specFactor',alpha_factor,'paramBounds',paramBounds);
                                            
% global eq                                            
[sim_info_eq, sim_struct_eq] = param_sweep_multi_v3([cw_index ir_index],functionPath,sweep_options{:},...
                                              'half_max_flag',true,...'cw',1,...
                                              'equilibrium_flag',true,'specFactor',alpha_factor,'paramBounds',paramBounds,'n_sim',20,'n_seeds',10);                                            

% s0-optimized
[sim_info_s0, sim_struct_s0] = param_sweep_multi_v3([cw_index sharp_right_norm_index],functionPath,sweep_options{:},...
                                              'half_max_flag',true,...'cw',1,...
                                              'equilibrium_flag',false,'specFactor',alpha_factor,'paramBounds',paramBounds,'n_sim',1);
                                            
% f0-optimized
[sim_info_f0, sim_struct_f0] = param_sweep_multi_v3([cw_index spec_index],functionPath,sweep_options{:},...
                                              'half_max_flag',true,...'cw',1,...
                                              'equilibrium_flag',false,'specFactor',alpha_factor,'paramBounds',paramBounds,'n_sim',1);   
                                            
% generic equilibrium                                            
% sweep vs rante entryopy, which should be largely uncorrelated with IR
[sim_info_eq0, sim_struct_eq0] = param_sweep_multi_v3([cw_index rate_ent_index],functionPath,sweep_options{:},...
                                              'half_max_flag',true,...'cw',1,...
                                              'equilibrium_flag',true,'specFactor',alpha_factor,'paramBounds',paramBounds);                                            
toc                                            
%% Identify optimal performers as a function of cw
cw_vec = logspace(-1,log10(alpha_factor^4),51);
n_eq_keep = 2e3;
n_neq_keep = 5e1; % per c value
rng(143);
% extract rate and metric arrays
metric_array_neq = vertcat(sim_struct_neq.metric_array);
cw_vec_neq = 10.^metric_array_neq(:,cw_index);
sharp_vec_neq = metric_array_neq(:,sharpness_index);
rate_array_neq = vertcat(sim_struct_neq.rate_array);

metric_array_eq = vertcat(sim_struct_eq.metric_array);
cw_vec_eq = 10.^metric_array_eq(:,cw_index);
sharp_vec_eq = metric_array_eq(:,sharpness_index);
rate_array_eq = vertcat(sim_struct_eq.rate_array);

metric_array_eq0 = vertcat(sim_struct_eq0.metric_array);
cw_vec_eq0 = 10.^metric_array_eq0(:,cw_index);
sharp_vec_eq0 = metric_array_eq0(:,sharpness_index);
rate_array_eq0 = vertcat(sim_struct_eq0.rate_array);
eq0_keep_indices = randsample(find(sharp_vec_eq0>=0),min([n_eq_keep,sum(sharp_vec_eq0>=0)]),false);

metric_array_s0 = vertcat(sim_struct_s0.metric_array);
cw_vec_s0 = 10.^metric_array_s0(:,cw_index);
sharp_vec_s0 = metric_array_s0(:,sharpness_index);
rate_array_s0 = vertcat(sim_struct_s0.rate_array);

metric_array_f0 = vertcat(sim_struct_f0.metric_array);
cw_vec_f0 = 10.^metric_array_f0(:,cw_index);
sharp_vec_f0 = metric_array_f0(:,sharpness_index);
rate_array_f0 = vertcat(sim_struct_f0.rate_array);


% initialize arrays
best_ir_indices_eq = NaN(1,length(cw_vec)-1);
best_ir_indices_neq = cell(1,length(cw_vec)-1);
best_s0_indices = NaN(1,length(cw_vec)-1);
best_f0_indices = NaN(1,length(cw_vec)-1);

best_ir_values_eq = NaN(1,length(cw_vec)-1);
best_ir_values_neq = cell(1,length(cw_vec)-1);
best_s0_values = NaN(1,length(cw_vec)-1);
best_f0_values = NaN(1,length(cw_vec)-1);
cw_values_neq = cell(1,length(cw_vec)-1);

for c = 1:length(cw_vec)-1
  
    % obtain filter vecs
    cw_filter_eq = cw_vec_eq >= cw_vec(c) & cw_vec_eq < cw_vec(c+1) & sharp_vec_eq >=0;
    cw_filter_neq = cw_vec_neq >= cw_vec(c) & cw_vec_neq < cw_vec(c+1) & sharp_vec_neq >=0;
    cw_filter_s0 = cw_vec_s0 >= cw_vec(c) & cw_vec_s0 < cw_vec(c+1) & sharp_vec_s0 >=0;
    cw_filter_f0 = cw_vec_f0 >= cw_vec(c) & cw_vec_f0 < cw_vec(c+1) & sharp_vec_f0 >=0;
  
    % find optima
    [best_ir_values_eq(c), best_ir_indices_eq(c)] = nanmax(metric_array_eq(:,ir_index).*cw_filter_eq);    
%     [best_ir_values_neq(c), best_ir_indices_neq(c)] = nanmax(metric_array_neq(:,ir_index).*cw_filter_neq);    
    [best_s0_values(c), best_s0_indices(c)] = nanmax(metric_array_s0(:,sharp_right_norm_index).*cw_filter_s0);
    [best_f0_values(c), best_f0_indices(c)] = nanmax(metric_array_f0(:,spec_index).*cw_filter_f0);
    
    % find best N
    ir_vec_neq = metric_array_neq(cw_filter_neq,ir_index);
    ir_best = nanmax(ir_vec_neq);
    options = find((metric_array_neq(:,ir_index)>=0.95*ir_best&cw_filter_neq));
    if length(options)>1
        best_ir_indices_neq{c} = randsample(options,min([n_neq_keep length(options)]),false);
    else
        best_ir_indices_neq{c} = options;
    end
    best_ir_values_neq{c} =  metric_array_neq(best_ir_indices_neq{c},ir_index);
    cw_values_neq{c} = cw_vec_neq(best_ir_indices_neq{c});% repelem(cw_vec(c+1),length(best_ir_values_neq{c}))';
end    
best_ir_indices_neq = vertcat(best_ir_indices_neq{:});
best_ir_values_neq = vertcat(best_ir_values_neq{:});
cw_vec_neq_plot = vertcat(cw_values_neq{:});

% add functions to path
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));

% set basic parameters
n_perturbations = 1e2;
m_factor = logspace(0,log10(alpha_factor),n_perturbations);
alpha_adjusted = alpha_factor ./ m_factor;

% generate arrays of optimal rates

% eq global
ir_rate_array_eq = repmat(rate_array_eq(best_ir_indices_eq,:), 1,1,length(m_factor));
ir_rate_array_eq(:,sim_info_eq.unbindingFlags==1,:) = ir_rate_array_eq(:,sim_info_eq.unbindingFlags==1,:) .* reshape(m_factor,1,1,[]);
ir_rate_array_eq(:,sim_info_eq.b_index,:) = repmat(reshape(alpha_adjusted,1,1,[]),size(ir_rate_array_eq,1),1,1);

% neq global
ir_rate_array_neq = repmat(rate_array_neq(best_ir_indices_neq,:),1,1, length(m_factor));
ir_rate_array_neq(:,sim_info_neq.unbindingFlags==1,:) = ir_rate_array_neq(:,sim_info_eq.unbindingFlags==1,:) .* reshape(m_factor,1,1,[]);
ir_rate_array_neq(:,sim_info_neq.b_index,:) = repmat(reshape(alpha_adjusted,1,1,[]),size(ir_rate_array_neq,1),1,1);

% s0
s0_rate_array = repmat(rate_array_s0(best_s0_indices,:),1,1, length(m_factor));
s0_rate_array(:,sim_info_neq.unbindingFlags==1,:) = s0_rate_array (:,sim_info_eq.unbindingFlags==1,:) .* reshape(m_factor,1,1,[]);
s0_rate_array(:,sim_info_neq.b_index,:) = repmat(reshape(alpha_adjusted,1,1,[]),size(s0_rate_array ,1),1,1);

% f0
f0_rate_array = repmat(rate_array_f0(best_f0_indices,:),1,1, length(m_factor));
f0_rate_array(:,sim_info_neq.unbindingFlags==1,:) = f0_rate_array(:,sim_info_eq.unbindingFlags==1,:) .* reshape(m_factor,1,1,[]);
f0_rate_array(:,sim_info_neq.b_index,:) = repmat(reshape(alpha_adjusted,1,1,[]),size(f0_rate_array,1),1,1);

% now do this for the equilibrium "hord"
eq_rate_array = repmat(rate_array_eq0(eq0_keep_indices,:), 1,1,length(m_factor));
eq_rate_array(:,sim_info_eq.unbindingFlags==1,:) = eq_rate_array(:,sim_info_eq.unbindingFlags==1,:) .* reshape(m_factor,1,1,[]);
eq_rate_array(:,sim_info_eq.b_index,:) = repmat(reshape(alpha_adjusted,1,1,[]),size(eq_rate_array,1),1,1);

%%
% calculate predicted observed sharpness and production for each kind of 
% network under different perturbation strengths

sharpness_vec_eq = NaN(size(ir_rate_array_eq,1),size(ir_rate_array_eq,3));
sharpness_vec_neq = NaN(size(ir_rate_array_neq,1),size(ir_rate_array_neq,3));
sharpness_vec_s0 = NaN(size(s0_rate_array,1),size(s0_rate_array,3));
sharpness_vec_f0 = NaN(size(f0_rate_array,1),size(f0_rate_array,3));
sharpness_vec_eq_full = NaN(size(eq_rate_array,1),size(eq_rate_array,3));

rate_vec_eq = NaN(size(ir_rate_array_eq,1),size(ir_rate_array_eq,3));
rate_vec_neq = NaN(size(ir_rate_array_neq,1),size(ir_rate_array_neq,3));
rate_vec_s0 = NaN(size(s0_rate_array,1),size(s0_rate_array,3));
rate_vec_f0 = NaN(size(f0_rate_array,1),size(f0_rate_array,3));
rate_vec_eq_full = NaN(size(eq_rate_array,1),size(eq_rate_array,3));

for p = 1:size(ir_rate_array_eq,3)
  
    % extract rates
    eq_rates = ir_rate_array_eq(:,:,p);    
    neq_rates = ir_rate_array_neq(:,:,p);   
    s0_rates = s0_rate_array(:,:,p);   
    f0_rates = f0_rate_array(:,:,p);   
    eq_rates_full = eq_rate_array(:,:,p);   
    
    % convert to cell array and run calculations
    paramCellEq = mat2cell(eq_rates,size(eq_rates,1),ones(1,size(eq_rates,2)));
    paramCellNeq = mat2cell(neq_rates,size(neq_rates,1),ones(1,size(neq_rates,2)));
%     paramCellS0 = mat2cell(s0_rates,size(s0_rates,1),ones(1,size(s0_rates,2)));
%     paramCellF0 = mat2cell(f0_rates,size(f0_rates,1),ones(1,size(f0_rates,2)));
    paramCellEqFull = mat2cell(eq_rates_full,size(eq_rates_full,1),ones(1,size(eq_rates_full,2)));
    
    % production rate 
    rate_vec_eq(:,p) = productionRateFunction(paramCellEq{:});
    rate_vec_neq(:,p) = productionRateFunction(paramCellNeq{:});
%     rate_vec_s0(:,p) = productionRateFunction(paramCellS0{:});
%     rate_vec_f0(:,p) = productionRateFunction(paramCellF0{:});
    rate_vec_eq_full(:,p) = productionRateFunction(paramCellEqFull{:});
    
    % sharpness
    sharpness_vec_eq(:,p) = sharpnessFunction(paramCellEq{:});
    sharpness_vec_neq(:,p) = sharpnessFunction(paramCellNeq{:});    
%     sharpness_vec_s0(:,p) = sharpnessFunction(paramCellS0{:});    
%     sharpness_vec_f0(:,p) = sharpnessFunction(paramCellF0{:});    
    sharpness_vec_eq_full(:,p) = sharpnessFunction(paramCellEqFull{:});    
   
end


%% Let's see what the predicted sharpness shift looks like for different
% designate index of perturbation strength to use for plots
m_ind = 16;
close all

% normalize arrays
sharpness_array_norm_eq = sharpness_vec_eq ./ sharpness_vec_eq(:,1);
sharpness_array_norm_neq = sharpness_vec_neq ./ sharpness_vec_neq(:,1);
sharpness_array_norm_s0 = sharpness_vec_s0 ./ sharpness_vec_s0(:,1);
sharpness_array_norm_f0 = sharpness_vec_f0 ./ sharpness_vec_f0(:,1);
sharpness_array_norm_eq_full = sharpness_vec_eq_full ./ sharpness_vec_eq_full(:,1);

rate_array_norm_eq = rate_vec_eq./ rate_vec_eq(:,1);
rate_array_norm_neq = rate_vec_neq ./ rate_vec_neq(:,1);
rate_array_norm_s0 = rate_vec_s0 ./ rate_vec_s0(:,1);
rate_array_norm_f0 = rate_vec_f0./ rate_vec_f0(:,1);
rate_array_norm_eq_full = rate_vec_eq_full./ rate_vec_eq_full(:,1);



cw_vec_plot = cw_vec(2:end);
sharp_shift_fig = figure('Position',[100 100 512 256]);
cmap = brewermap([],'Set2');
hold on
scatter(cw_vec_eq0(eq0_keep_indices),sharpness_array_norm_eq_full(:,m_ind).*m_factor(m_ind),'s','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k')
scatter(cw_vec_neq_plot,sharpness_array_norm_neq(:,m_ind).*m_factor(m_ind),40,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k')
% scatter(cw_vec_plot,sharpness_array_norm_eq(:,m_ind),markerSize,'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k')
% scatter(cw_vec_plot,sharpness_array_norm_s0(:,m_ind),markerSize,'MarkerFaceColor',cmap(5,:),'MarkerEdgeColor','k')
% scatter(cw_vec_plot,sharpness_array_norm_f0(:,m_ind),markerSize,'MarkerFaceColor',cmap(4,:),'MarkerEdgeColor','k')
% set(gca,'yscale','log')
set(gca,'xscale','log')
xlabel('relative non-cognate factor concentration (c_w / c_r)');
ylabel('normalized sharpness shift (s^*/s \times k_-^*/k_-^0)')
xlim([1e0 alpha_factor^4])
ylim([.2 2])

set(gca,'xtick',[1 alpha_factor^1 alpha_factor^2 alpha_factor^3 alpha_factor^4]);%,'xticklabels',{'\alpha^0','\alpha^1','\alpha^2','\alpha^3','\alpha^4'})

set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
sharp_shift_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
% saveas(sharp_shift_fig,[FigPath 'sharp_fold_vs_cw.png'])
% saveas(sharp_shift_fig,[FigPath 'sharp_fold_vs_cw.pdf'])

% rate shift
rate_shift_fig = figure('Position',[100 100 512 256]);
cmap = brewermap([],'Set2');
hold on
scatter(cw_vec_eq0(eq0_keep_indices),rate_array_norm_eq_full(:,m_ind),'s','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k')
scatter(cw_vec_neq_plot,rate_array_norm_neq(:,m_ind),40,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3)
% scatter(cw_vec_plot,rate_array_norm_s0(:,m_ind),markerSize,'MarkerFaceColor',cmap(5,:),'MarkerEdgeColor','k')
% scatter(cw_vec_plot,rate_array_norm_f0(:,m_ind),markerSize,'MarkerFaceColor',cmap(4,:),'MarkerEdgeColor','k')
% set(gca,'yscale','log')
set(gca,'xscale','log')
xlabel('relative non-cognate factor concentration (c_w / c_r)');
ylabel('rate change (r^*/r)')

set(gca,'xtick',[1 alpha_factor^1 alpha_factor^2 alpha_factor^3 alpha_factor^4]);%,'xticklabels',{'\alpha^0','\alpha^1','\alpha^2','\alpha^3','\alpha^4'})
% legend([p2 p1 s1],'sharpness optimized','fidelity optimized','global optimum')
xlim([1e0 alpha_factor^4])
set(gca,'FontSize',14)
ylim([.3 1.1])
set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
set(gca,'yscale','log')
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
        scatter(rate_array_norm_eq_full(cw_eq0_indices==c,m_ind_2),sharpness_array_norm_eq_full(cw_eq0_indices==c,m_ind_2).*m_factor(m_ind_2),[],...
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
