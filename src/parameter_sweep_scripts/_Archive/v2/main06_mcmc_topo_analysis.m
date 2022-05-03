% script to analyze topology of optimal transcription networksnetworks
clear
close all
% specify paths
addpath('../utilities')
readPath = '../../out/parameter_search/';
mkdir(readPath)
% specify model hyperparameters
equilibrium_flag = 0; % constrained to equilibrium?
metric_names = {'cooperativity', 'information_rate', 'critical_time',...
    'production_sign', 'production_rate'};
metric_indices = [1,2];
% set save name
mcmc_master_struct = struct;
for i = 1:numel(metric_indices)
    read_name = ['mcmc_opt_results_eq' num2str(equilibrium_flag) '_' metric_names{metric_indices(1)} '.mat'];
    load([readPath read_name])
    mcmc_master_struct(i).mcmc_results = mcmc_results;
end
% now load parameter sweep results
read_name = ['edge_sim_results_eq' num2str(equilibrium_flag) '_' ...
    metric_names{metric_indices(1)} '_' metric_names{metric_indices(2)} '.mat'];
load([readPath read_name])
% define figure path
figPath = ['../../fig/topology_scatters_eq' num2str(equilibrium_flag) '/'];
mkdir(figPath)

% define color map
cmap = brewermap(128,'RdYlBu');

% iterate through archetypes and analyze topologies
for i = 1:numel(mcmc_results)
    % extract param sweep results
    rate_array_sweep = simulation_results(i).rate_array;
    pd_rate_vec_sweep = mcmc_results(i).metric_array.production_rate;
    cycle_flux_vec_sweep = mcmc_results(i).cycle_flux_vec;
    % extract fidelity MCMC results
    rate_array_fidelity = mcmc_master_struct(1).mcmc_results(i).rate_array;
    pd_rate_vec_fidelity = mcmc_master_struct(1).mcmc_results(i).metric_array.production_rate;
    cycle_flux_vec_fidelity = mcmc_master_struct(1).mcmc_results(i).cycle_flux_vec;
    fidelity_vec = mcmc_master_struct(1).mcmc_results(i).metric_array.cooperativity;
    fid_real = mcmc_master_struct(2).mcmc_results(i).real_flag_vec;
    % filter for top 10% performing models
    thresh_fid = prctile(fidelity_vec,95);
    fidelity_ids = find(fidelity_vec>thresh_fid & fid_real);
    % extract information rate MCMC results
    rate_array_information = mcmc_master_struct(2).mcmc_results(i).rate_array;
    pd_rate_vec_information = mcmc_master_struct(2).mcmc_results(i).metric_array.production_rate;
    cycle_flux_vec_information = mcmc_master_struct(2).mcmc_results(i).cycle_flux_vec;
    information_vec = mcmc_master_struct(2).mcmc_results(i).metric_array.information_rate;            
    info_real = mcmc_master_struct(2).mcmc_results(i).real_flag_vec;
    % filter for top 10% performing models
    thresh_info = prctile(information_vec,95);
    info_ids = find(information_vec>thresh_info & info_real);
    % get PCA for sweep results       
    [coeff, score_mat] = pca(rate_array_sweep);  
    % project other rate arrays across coefficients
    score_fidelity = (rate_array_fidelity-nanmean(rate_array_fidelity)) * inv(coeff');
    score_information = (rate_array_information-nanmean(rate_array_information)) * inv(coeff');
%     score_sweep = (rate_array_sweep-nanmean(rate_array_sweep)) * inv(coeff');
    % visualize first two components
    pca_scatter = figure;
    hold on
    scatter(score_mat(:,1),score_mat(:,2),20,'MarkerFaceColor',[.6 .6 .6],...
        'MarkerFaceAlpha',.05,'MarkerEdgeAlpha',0)
    scatter(score_fidelity(fidelity_ids,1),score_fidelity(fidelity_ids,2),20,'MarkerFaceColor',cmap(1,:),...
        'MarkerFaceAlpha',.05,'MarkerEdgeAlpha',0)
    scatter(score_information(info_ids,1),score_information(info_ids,2),20,'MarkerFaceColor',cmap(end,:),...
        'MarkerFaceAlpha',.05,'MarkerEdgeAlpha',0)
    xlabel('1st principal component')
    ylabel('2nd principal component')
    grid on
%     legend('all points','optimal points')
    saveas(pca_scatter,[figPath 'md' num2str(i) '_pca_scatter.tif'])
       
    cycle_pd_fig = figure;
    hold on
    scatter(abs(cycle_flux_vec_sweep),pd_rate_vec_sweep,20,'MarkerFaceColor',[.6 .6 .6],...
        'MarkerFaceAlpha',.05,'MarkerEdgeAlpha',0)
    scatter(abs(cycle_flux_vec_fidelity(fidelity_ids)),pd_rate_vec_fidelity(fidelity_ids)...
        ,20,'MarkerFaceColor',cmap(1,:),'MarkerFaceAlpha',.05,'MarkerEdgeAlpha',0)
    scatter(abs(cycle_flux_vec_information(info_ids)),pd_rate_vec_information(info_ids)...
        ,20,'MarkerFaceColor',cmap(end,:),'MarkerFaceAlpha',.05,'MarkerEdgeAlpha',0)
    set(gca,'YScale','log','XScale','log')
    xlabel('cycle flux magnitude (s^{-1})')
    ylabel('mRNA production rate')
    grid on
	saveas(cycle_pd_fig,[figPath 'md' num2str(i) '_cycle_production.tif'])
end
    
    



