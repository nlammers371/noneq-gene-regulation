% Plot results for sharpness vs precision for higher order models
clear 
close all
addpath(genpath('../utilities/'))

% %%%%%%%%%%%%%%%%  set relative read and write paths %%%%%%%%%%%%%%%%%%%%
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\';
DataPath = [DropboxFolder  'SweepOutput\sweeps05B_info_vs_cw' filesep ];
DataPathHM = [DropboxFolder  'SweepOutput\sweeps05_info_vs_cw_hm' filesep ];
FigPath = [DropboxFolder '\manuscript\experimental_signatures' filesep];
mkdir(FigPath);

% %%%%%%%%%%%%%%%%  Set plot parameters and constants %%%%%%%%%%%%%%%%%%%%
n_plot = 3e3; % number of points to plot
n_samp = 5e4; % number of points to plot
markerAlpha = 0.5; % marker transparency
markerSize = 75; % marker size
bit_factor = log2(exp(1));
alpha_factor = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%% Compare "equilibrium" strategy of adding binding sites with 
%%% "Non-equilibrium" strategy of adding locus conformations

% get metric names for numeric sweeps
[~,~,metric_names_num] = calculateMetricsNumeric_v3([]);

ir_index = find(strcmp(metric_names_num,'IR'));
cw_index = find(strcmp(metric_names_num,'CW'));
rate_index = find(strcmp(metric_names_num,'ProductionRate'));
sharpness_index = find(strcmp(metric_names_num,'Sharpness'));


%%%%%%%%%%%%%%%%%%%%%%%%5
%%% Load multi bs results
% get list of sweep results files with only 1 genera TF reaction
multi_bs_sweep_files_eq = dir([DataPath 'sweep_results*eq*']);
multi_bs_info_files_eq = dir([DataPath 'sweep_info*eq*']);

multi_bs_sweep_files_neq = dir([DataPath 'sweep_results*_cw1.mat']);
multi_bs_info_files_neq = dir([DataPath 'sweep_info*cw1.mat']);

% load
master_struct_multi = struct;
for f = 1:length(multi_bs_sweep_files_eq)
  
    % load eq files
    load([DataPath multi_bs_sweep_files_eq(f).name])
    load([DataPath multi_bs_info_files_eq(f).name])
    
    master_struct_multi(f).sweep_results_eq = sim_results;
    master_struct_multi(f).sweep_info_eq = sim_info;
        
    % load neq files
    load([DataPath multi_bs_sweep_files_neq(f).name])
    load([DataPath multi_bs_info_files_neq(f).name])
    
    master_struct_multi(f).sweep_results = sim_results;
    master_struct_multi(f).sweep_info = sim_info;
end    

%% Identify optimal performers as a function of cw
cw_bins = logspace(0,5,51);
n_eq_keep = 5e4;
n_keep = 1e3; % per c value
rng(143);
r1 = 0.9;
r0 = 0.1;

n_perturbations = 1e2;
% m_factor = logspace(0,log10(alpha_factor),n_perturbations);
m_factor = [1 1.25 2];
alpha_adjusted = alpha_factor ./ m_factor;

% indicate indices of models to use
g_index = 1;
bs_index = 2;

% extract rate and metric arrays
metric_array_neq_g = vertcat(master_struct_multi(g_index).sweep_results.metric_array);
cw_vec_neq_g = 10.^metric_array_neq_g(:,cw_index);
ir_vec_neq_g = metric_array_neq_g(:,ir_index);
sharp_vec_neq_g = metric_array_neq_g(:,sharpness_index);
rate_vec_neq_g = metric_array_neq_g(:,rate_index);
rate_array_neq_g = vertcat(master_struct_multi(g_index).sweep_results.rate_array);

metric_array_eq_g = vertcat(master_struct_multi(g_index).sweep_results_eq.metric_array);
cw_vec_eq_g = 10.^metric_array_eq_g(:,cw_index);
sharp_vec_eq_g = metric_array_eq_g(:,sharpness_index);
ir_vec_eq_g = metric_array_eq_g(:,ir_index);
rate_vec_eq_g = metric_array_eq_g(:,rate_index);
rate_array_eq_g = vertcat(master_struct_multi(g_index).sweep_results_eq.rate_array);

metric_array_neq_bs = vertcat(master_struct_multi(bs_index).sweep_results.metric_array);
cw_vec_neq_bs = 10.^metric_array_neq_bs(:,cw_index);
sharp_vec_neq_bs = metric_array_neq_bs(:,sharpness_index);
ir_vec_neq_bs = metric_array_neq_bs(:,ir_index);
rate_vec_neq_bs = metric_array_neq_bs(:,rate_index);
rate_array_neq_bs = vertcat(master_struct_multi(bs_index).sweep_results.rate_array);

metric_array_eq_bs = vertcat(master_struct_multi(bs_index).sweep_results_eq.metric_array);
cw_vec_eq_bs = 10.^metric_array_eq_bs(:,cw_index);
sharp_vec_eq_bs = metric_array_eq_bs(:,sharpness_index);
ir_vec_eq_bs = metric_array_eq_bs(:,ir_index);
rate_vec_eq_bs = metric_array_eq_bs(:,rate_index);
rate_array_eq_bs = vertcat(master_struct_multi(bs_index).sweep_results_eq.rate_array);

% initialize arrays
best_ir_index_cell_g_neq = cell(1,length(cw_bins)-1);
best_ir_value_cell_g_neq = cell(1,length(cw_bins)-1);
best_ir_index_cell_g_eq = cell(1,length(cw_bins)-1);
best_ir_value_cell_g_eq = cell(1,length(cw_bins)-1);

best_ir_index_cell_bs_neq = cell(1,length(cw_bins)-1);
best_ir_value_cell_bs_neq = cell(1,length(cw_bins)-1);
best_ir_index_cell_bs_eq = cell(1,length(cw_bins)-1);
best_ir_value_cell_bs_eq = cell(1,length(cw_bins)-1);

for c = 1:length(cw_bins)-1  
    % obtain filter vecs 
    cw_filter_g_neq = cw_vec_neq_g >= cw_bins(c) & cw_vec_neq_g < cw_bins(c+1) & sharp_vec_neq_g >=0;
    rate_filter_g_neq = rate_vec_neq_g >= r0 & rate_vec_neq_g <=r1 & ir_vec_neq_g>=1e-5;
    
    cw_filter_g_eq = cw_vec_eq_g >= cw_bins(c) & cw_vec_eq_g < cw_bins(c+1) & sharp_vec_eq_g >=0;
    rate_filter_g_eq = rate_vec_eq_g >= r0 & rate_vec_eq_g <=r1 & ir_vec_eq_g>=1e-5;
    
    cw_filter_bs_neq = cw_vec_neq_bs >= cw_bins(c) & cw_vec_neq_bs < cw_bins(c+1) & sharp_vec_neq_bs >=0;
    rate_filter_bs_neq = rate_vec_neq_bs >= r0 & rate_vec_neq_bs <=r1 & ir_vec_neq_bs>=1e-5;
    
    cw_filter_bs_eq = cw_vec_eq_bs >= cw_bins(c) & cw_vec_eq_bs < cw_bins(c+1) & sharp_vec_eq_bs >=0;
    rate_filter_bs_eq = rate_vec_eq_bs >= r0 & rate_vec_eq_bs <=r1& ir_vec_eq_bs>=1e-8;
    
    % assign percentile scores to remaining
        
    keep_indices = find(cw_filter_g_neq&rate_filter_g_neq);
    ir_max_neq_g = nanmax(metric_array_neq_g(keep_indices,ir_index));
    best_ir_index_cell_g_neq{c} = keep_indices;
    best_ir_value_cell_g_neq{c} = metric_array_neq_g(keep_indices,ir_index)/ir_max_neq_g;
    
    keep_indices = find(cw_filter_g_eq&rate_filter_g_eq);
    best_ir_index_cell_g_eq{c} = keep_indices;
    best_ir_value_cell_g_eq{c} = metric_array_eq_g(keep_indices,ir_index)/ir_max_neq_g;
    
    keep_indices = find(cw_filter_bs_neq&rate_filter_bs_neq);
    ir_max_neq_bs = nanmax(metric_array_neq_bs(keep_indices,ir_index)); 
    best_ir_index_cell_bs_neq{c} = keep_indices;
    best_ir_value_cell_bs_neq{c} = metric_array_neq_bs(keep_indices,ir_index)/ir_max_neq_bs;
    
    keep_indices = find(cw_filter_bs_eq&rate_filter_bs_eq);
    best_ir_index_cell_bs_eq{c} = keep_indices;
    best_ir_value_cell_bs_eq{c} = metric_array_eq_bs(keep_indices,ir_index)/ir_max_neq_bs;
end    

% % add functions to path
% rmpath(genpath('../utilities/metricFunctions/'));
% addpath(genpath(functionPath));

% set basic parameters
sim_info_g_neq = master_struct_multi(g_index).sweep_info;
sim_info_g_neq.unbindingFlags = strcmp(sim_info_g_neq.sweepVarStrings,'km');
sim_info_g_eq = master_struct_multi(1).sweep_info;
sim_info_g_eq.unbindingFlags = strcmp(sim_info_g_eq.sweepVarStrings,'km');
sim_info_bs_neq = master_struct_multi(bs_index).sweep_info;
sim_info_bs_neq.unbindingFlags = strcmp(sim_info_bs_neq.sweepVarStrings,'km');
sim_info_bs_neq.coopFlags = strcmp(sim_info_bs_neq.sweepVarStrings,'wmp');


% generate arrays of optimal rates

% neq multi LC
ir_keep_indices_g_neq = vertcat(best_ir_index_cell_g_neq{:});
ir_keep_values_g_neq = vertcat(best_ir_value_cell_g_neq{:});
ir_rate_array_g_neq = repmat(rate_array_neq_g(ir_keep_indices_g_neq,:),1,1, length(m_factor));
% ir_rate_array_g_neq(:,sim_info_g_neq.unbindingFlags==1,:) = ir_rate_array_g_neq(:,sim_info_g_neq.unbindingFlags==1,:) .* reshape(m_factor,1,1,[]);
% ir_rate_array_g_neq(:,sim_info_g_neq.a_index,:) = repmat(reshape(alpha_adjusted,1,1,[]),size(ir_rate_array_g_neq,1),1,1);

% eq multi LC
ir_keep_indices_g_eq = vertcat(best_ir_index_cell_g_eq{:});
ir_keep_values_g_eq = vertcat(best_ir_value_cell_g_eq{:});
ir_rate_array_g_eq = repmat(rate_array_eq_g(ir_keep_indices_g_eq,:),1,1, length(m_factor));
% ir_rate_array_g_eq(:,sim_info_g_eq.unbindingFlags==1,:) = ir_rate_array_g_eq(:,sim_info_g_eq.unbindingFlags==1,:) .* reshape(m_factor,1,1,[]);
% ir_rate_array_g_eq(:,sim_info_g_eq.a_index,:) = repmat(reshape(alpha_adjusted,1,1,[]),size(ir_rate_array_g_eq,1),1,1);

% neq multi BS
ir_keep_indices_bs_neq = vertcat(best_ir_index_cell_bs_neq{:});
ir_keep_values_bs_neq = vertcat(best_ir_value_cell_bs_neq{:});
ir_rate_array_bs_neq = repmat(rate_array_neq_bs(ir_keep_indices_bs_neq,:),1,1, length(m_factor));
% ir_rate_array_bs_neq(:,sim_info_bs_neq.unbindingFlags==1,:) = ir_rate_array_bs_neq(:,sim_info_bs_neq.unbindingFlags==1,:) .* reshape(m_factor,1,1,[]);
% ir_rate_array_bs_neq(:,sim_info_bs_neq.coopFlags==1,:) = 1;
% ir_rate_array_bs_neq(:,sim_info_bs_neq.a_index,:) = repmat(reshape(alpha_adjusted,1,1,[]),size(ir_rate_array_bs_neq,1),1,1);

% eq multi BS
ir_keep_indices_bs_eq = vertcat(best_ir_index_cell_bs_eq{:});
ir_keep_values_bs_eq = vertcat(best_ir_value_cell_bs_eq{:});
ir_rate_array_bs_eq = repmat(rate_array_eq_bs(ir_keep_indices_bs_eq,:),1,1, length(m_factor));
% ir_rate_array_bs_eq(:,sim_info_bs_neq.unbindingFlags==1,:) = ir_rate_array_bs_eq(:,sim_info_bs_neq.unbindingFlags==1,:) .* reshape(m_factor,1,1,[]);
% ir_rate_array_bs_eq(:,sim_info_bs_neq.coopFlags==1,:) = 1;
% ir_rate_array_bs_eq(:,sim_info_bs_neq.a_index,:) = repmat(reshape(alpha_adjusted,1,1,[]),size(ir_rate_array_bs_eq,1),1,1);

%% calculate predicted observed sharpness and production for each kind of 

% state probabilities
numerical_precision = 5;
    
if true%~exist([FigPath 'mutation_results.mat'])
    % network under different perturbation strengths
    sharpness_array_g_neq = NaN(size(ir_rate_array_g_neq,1),size(ir_rate_array_g_neq,3));
    sharpness_array_g_eq = NaN(size(ir_rate_array_g_eq,1),size(ir_rate_array_g_eq,3));
    sharpness_array_bs_neq = NaN(size(ir_rate_array_bs_neq,1),size(ir_rate_array_bs_neq,3));
    sharpness_array_bs_eq = NaN(size(ir_rate_array_bs_eq,1),size(ir_rate_array_bs_eq,3));

    rate_array_g_neq = NaN(size(ir_rate_array_g_neq,1),size(ir_rate_array_g_neq,3));
    rate_array_g_eq = NaN(size(ir_rate_array_g_eq,1),size(ir_rate_array_g_eq,3));
    rate_array_bs_neq = NaN(size(ir_rate_array_bs_neq,1),size(ir_rate_array_bs_neq,3));
    rate_array_bs_eq = NaN(size(ir_rate_array_bs_eq,1),size(ir_rate_array_bs_eq,3));

    for p = 1:size(ir_rate_array_g_neq,3)
        tic

        % extract rates   
        neq_rates_g = ir_rate_array_g_neq(:,:,p);
        eq_rates_g = ir_rate_array_g_eq(:,:,p);
        neq_rates_bs = ir_rate_array_bs_neq(:,:,p);   
        eq_rates_bs = ir_rate_array_bs_eq(:,:,p);       

        %%%%%%%%%%%%%%%%%%%
        % NEQ G
        for i = 3:4
            if i == 1
                % set function path
                sweepInfo = master_struct_multi(g_index).sweep_info;             
                rates = neq_rates_g;
            elseif i == 2
                sweepInfo = master_struct_multi(g_index).sweep_info_eq;             
                rates = eq_rates_g;
            elseif i == 3
                sweepInfo = master_struct_multi(bs_index).sweep_info;             
                rates = neq_rates_bs;
            else
                sweepInfo = master_struct_multi(bs_index).sweep_info_eq;             
                rates = eq_rates_bs;
            end
            %
            % extract full transition matrix
            RSymFull = sweepInfo.RSymFull;            
            nStates = size(RSymFull,1);
            
            % find entries to alter
            [ub_cols, si] = sort([2:6:nStates 3:6:nStates]);
            ub_rows = [1:6:nStates 4:6:nStates];
            ub_rows = ub_rows(si);
            lin_ids = sub2ind(size(RSymFull),ub_rows,ub_cols);
            
            % alter entries
            RSymFullMut = RSymFull;
            RSymFullMut(lin_ids) = RSymFullMut(lin_ids)*m_factor(p);
            
            % update diagonal terms
            ind_mat = eye(nStates);
            RSymFullMut(ind_mat==1) = 0;
            RSymFullMut(ind_mat==1) = -sum(RSymFullMut);            
            u_bound_vec = zeros(1, nStates);
            u_bound_vec(ub_cols) = 1;
            
            % condense as much as possible            
            full_state_array = [sweepInfo.n_g_engaged_full' sweepInfo.n_right_bound_full' sweepInfo.n_wrong_bound_full' u_bound_vec'];
            
            [unique_state_array,ia,~] = unique(full_state_array,'rows');
            
            RSymSmall = sym(zeros(size(unique_state_array,1)));
            for j = 1:size(unique_state_array,1)        
                state_ids = ismember(full_state_array,unique_state_array(j,:),'rows');%n_g_engaged==unique_state_array(i,1) & n_right_bound==unique_state_array(i,2) & n_wrong_bound==unique_state_array(i,3);
                to_col = sum(RSymFullMut(state_ids,:),1);
                RSymSmall(j,:) = to_col(ia);
            end
            
            % make matlab function
            tic
            RFunTemp = matlabFunction(RSymSmall,'vars',sweepInfo.sweepVarList,'Optimize',false);
            toc
            
            % generate activity vector
            n_g_engaged_small = sweepInfo.n_g_engaged_full(ia);
            active_state_filter = n_g_engaged_small==max(n_g_engaged_small);
            
            
%             functionPath = sweepInfo.functionPath; 
%             slashes = strfind(functionPath,'\');
%             simName = functionPath(slashes(end-1)+1:slashes(end)-1);
%             rmpath(genpath('../utilities/metricFunctions/'));
%             addpath(genpath(['../utilities/metricFunctions/numeric/' simName]));

            pd_temp = NaN(1,size(rates,1));
            s_temp = NaN(1,size(rates,1));
            for r = 1:size(rates,1)

                paramCell = mat2cell(rates(r,:),size(rates(r,:),1),ones(1,size(rates(r,:),2)));
                Q_num = RFunTemp(paramCell{:});

                ss_full = calculate_ss_num(Q_num,numerical_precision);  

                % get production rate
                pd_temp(r) = sum(ss_full(active_state_filter));

                % get sharpness            
                paramsC0 = rates(r,:);
                paramsC0(sweepInfo.cr_index) = sweepInfo.cr0;
                paramsC1 = rates(r,:);
                paramsC1(sweepInfo.cr_index) = sweepInfo.cr1;               
                paramCellC0 = mat2cell(paramsC0,size(paramsC0,1),ones(1,size(paramsC0,2)));
                paramCellC1 = mat2cell(paramsC1,size(paramsC1,1),ones(1,size(paramsC1,2)));

                Q_num_c0 = RFunTemp(paramCellC0{:});
                Q_num_c1 = RFunTemp(paramCellC1{:});

                % calculate probabilities                            
                state_probs_c0 = calculate_ss_num(Q_num_c0,numerical_precision)';
                state_probs_c1 = calculate_ss_num(Q_num_c1,numerical_precision)';

                ProductionRate1 = sum(state_probs_c1(:,active_state_filter),2);
                ProductionRate0 = sum(state_probs_c0(:,active_state_filter),2);   

                s_temp(r) = (ProductionRate1-ProductionRate0)./(sweepInfo.cr1-sweepInfo.cr0); 
            end

            % save output
            if i == 1
                sharpness_array_g_neq(:,p) = s_temp;
                rate_array_g_neq(:,p) = pd_temp;
            elseif i == 2
                sharpness_array_g_eq(:,p) = s_temp;
                rate_array_g_eq(:,p) = pd_temp;
            elseif i == 3
                sharpness_array_bs_neq(:,p) = s_temp;
                rate_array_bs_neq(:,p) = pd_temp;
            else           
                sharpness_array_bs_eq(:,p) = s_temp;
                rate_array_bs_eq(:,p) = pd_temp;
            end
        end
        toc
    end

    % save
    master_struct_multi(g_index).sharpness_array_neq = sharpness_array_g_neq;
    master_struct_multi(g_index).rate_array_neq = rate_array_g_neq;
    master_struct_multi(g_index).sharpness_array_eq = sharpness_array_g_eq;
    master_struct_multi(g_index).rate_array_eq = rate_array_g_eq;
    
    master_struct_multi(bs_index).sharpness_array_neq = sharpness_array_bs_neq;
    master_struct_multi(bs_index).rate_array_neq = rate_array_bs_neq;
    master_struct_multi(bs_index).sharpness_array_eq = sharpness_array_bs_eq;
    master_struct_multi(bs_index).rate_array_eq = rate_array_bs_eq;

    save([FigPath 'mutation_results.mat'],'master_struct_multi')
else
    sharpness_array_g_neq = master_struct_multi(g_index).sharpness_array_neq;
    rate_array_g_neq = master_struct_multi(g_index).rate_array_neq;
    sharpness_array_g_eq = master_struct_multi(g_index).sharpness_array_eq;
    rate_array_g_eq = master_struct_multi(g_index).rate_array_eq;
    
    sharpness_array_bs_neq = master_struct_multi(bs_index).sharpness_array_neq;
    rate_array_bs_neq = master_struct_multi(bs_index).rate_array_neq;
    sharpness_array_bs_eq = master_struct_multi(bs_index).sharpness_array_eq;
    rate_array_bs_eq = master_struct_multi(bs_index).rate_array_eq;
end
%%
% Let's see what the predicted sharpness shift looks like for different
% designate index of perturbation strength to use for plots
m_ind = 3;
close all

% normalize arrays
sharpness_array_norm_g_neq = sharpness_array_g_neq ./ sharpness_array_g_neq(:,1);
sharpness_array_norm_g_eq = sharpness_array_g_eq ./ sharpness_array_g_eq(:,1);
sharpness_array_norm_bs_neq = sharpness_array_bs_neq ./ sharpness_array_bs_neq(:,1);
sharpness_array_norm_bs_eq = sharpness_array_bs_eq ./ sharpness_array_bs_eq(:,1);

rate_array_norm_g_neq = rate_array_g_neq ./ rate_array_g_neq(:,1);
rate_array_norm_g_eq = rate_array_g_eq ./ rate_array_g_eq(:,1);
rate_array_norm_bs_neq = rate_array_bs_neq ./ rate_array_bs_neq(:,1);
rate_array_norm_bs_eq = rate_array_bs_eq ./ rate_array_bs_eq(:,1);

% sharpness_array_norm_g_neq = sharpness_array_g_neq - sharpness_array_g_neq(:,1);
% sharpness_array_norm_g_eq = sharpness_array_g_eq - sharpness_array_g_eq(:,1);
% sharpness_array_norm_bs_neq = sharpness_array_bs_neq - sharpness_array_bs_neq(:,1);
% sharpness_array_norm_bs_eq = sharpness_array_bs_eq - sharpness_array_bs_eq(:,1);
% 
% rate_array_norm_g_neq = rate_array_g_neq - rate_array_g_neq(:,1);

% rate_array_norm_g_eq = rate_array_g_eq - rate_array_g_eq(:,1);
% rate_array_norm_bs_neq = rate_array_bs_neq - rate_array_bs_neq(:,1);
% rate_array_norm_bs_eq = rate_array_bs_eq - rate_array_bs_eq(:,1);

cw_vec_plot = cw_bins(2:end);

%%
close all

sharp_shift_fig = figure;%('Position',[100 100 512 256]);
cmap = brewermap([],'Set2');
cmap2 = flipud(brewermap([],'Spectral'));
colormap(cmap2)
hold on
% scatter(cw_vec_neq_g(ir_keep_indices_g_neq),rate_array_norm_g_neq(:,m_ind),[],ir_keep_values_g_neq,'filled','o','MarkerEdgeColor','k','MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0)
% scatter(cw_vec_eq_g(ir_keep_indices_g_eq),rate_array_norm_g_eq(:,m_ind),'s','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k')
scatter(cw_vec_neq_bs(ir_keep_indices_bs_neq),rate_array_norm_bs_neq(:,m_ind),[],ir_keep_values_bs_neq,'filled','o','MarkerEdgeColor','k','MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0)
scatter(cw_vec_eq_bs(ir_keep_indices_bs_eq),rate_array_norm_bs_eq(:,m_ind),'s','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k')

set(gca,'xscale','log')
xlabel('relative non-cognate factor concentration (c / c)');
ylabel('normalized sharpness shift (s/s \times k/k)')
xlim([1e0 1e5])
% ylim([0 4])

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
% saveas(sharp_shift_fig,[FigPath 'sharp_fold_vs_cw.png'])
% saveas(sharp_shift_fig,[FigPath 'sharp_fold_vs_cw.pdf'])
%%
c_vec_neq = cw_vec_neq_bs(ir_keep_indices_bs_neq);
rate_vec_neq = rate_vec_neq_bs(ir_keep_indices_bs_neq);
b_vec_neq = rate_vec_neq.*(1-rate_vec_neq);

c_vec_eq = cw_vec_eq_bs(ir_keep_indices_bs_eq);
rate_vec_eq = rate_vec_eq_bs(ir_keep_indices_bs_eq);
b_vec_eq = rate_vec_eq.*(1-rate_vec_eq);

r_ft_neq = rate_vec_neq>=0&rate_vec_neq<=1;
r_ft_eq = rate_vec_eq>=0&rate_vec_eq<=1;

cval = 1e2;
figure;
hold on
% scatter(sharpness_array_bs_neq(r_ft_neq&c_vec_neq>=cval,1)./b_vec_neq(r_ft_neq&c_vec_neq>=cval),rate_array_norm_bs_neq(r_ft_neq&c_vec_neq>=cval,m_ind)*2,[],ir_keep_values_bs_neq(r_ft_neq&c_vec_neq>=cval)>=0.95,'filled','o','MarkerEdgeColor','k','MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0)
% scatter(sharpness_array_bs_eq(r_ft_eq&c_vec_eq>=cval,1)./b_vec_eq(r_ft_eq&c_vec_eq>=cval),rate_array_norm_bs_eq(r_ft_eq&c_vec_eq>=cval,m_ind)*2,'s','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k')
scatter(rate_array_norm_bs_neq(r_ft_neq&c_vec_neq>=cval,m_ind),sharpness_array_norm_bs_neq(r_ft_neq&c_vec_neq>=cval,m_ind)*1.25,[],ir_keep_values_bs_neq(r_ft_neq&c_vec_neq>=cval)>=0.9,'filled','o','MarkerEdgeColor','k','MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0)
scatter(rate_array_norm_bs_eq(r_ft_eq&c_vec_eq>=cval,m_ind),sharpness_array_norm_bs_eq(r_ft_eq&c_vec_eq>=cval,m_ind)*1.25,'s','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k')
% ylim([0 4])

%%
close all
% rate shift
rate_shift_fig = figure;%('Position',[100 100 512 256]);
cmap = brewermap([],'Set2');
cmap2 = flipud(brewermap([],'Spectral'));
colormap(cmap2)
ir_ft = ir_keep_values_g_neq>0.99;
hold on
scatter(cw_vec_neq_g(ir_keep_indices_g_neq),rate_array_norm_neq_full(:,m_ind),[],ir_keep_values_g_neq,'filled','o','MarkerEdgeColor','k','MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0)
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
cw_eq0_indices = discretize(cw_vec_eq0_plot,cw_bins);

phase_space_fig = figure;
hold on
s = [];
for m = 1:length(m_ind_vec)
    m_ind_2 = m_ind_vec(m);
    cmap = c_cell{m};
    r_vec = imgaussfilt(rate_array_norm_neq(:,m_ind_2),1);
    s_vec = imgaussfilt(sharpness_array_norm_g_neq(:,m_ind_2),1);
    for c = 1:length(cw_bins)-1
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
