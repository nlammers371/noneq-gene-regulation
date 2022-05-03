clear
close all
% specify paths
addpath('../utilities')
readPath = '../../out/parameter_search/';
figPath = '../../fig/parameter_search/';
suffix = 'cooperativity_information_rate';
figPath = [figPath suffix '/'];

% inititalize structure to save eq and non-eq results
master_struct = struct;
for i = 1:2
    load([readPath 'edge_sim_results_eq' num2str(i-1) '_' suffix '.mat'])
    master_struct(i).simulation_results = simulation_results;
end
% define color map
cmap = flipud(brewermap(128,'RdYlBu'));
% perform simple PCA analysis on non-eq solutions
% iterate through each network topology
for i = numel(simulation_results)-1
    % extract key quantities
    non_eq_rates = master_struct(1).simulation_results(i).rate_array;
    information_rate = master_struct(1).simulation_results(i).metric_array(:,2);
    fidelity_vec = master_struct(1).simulation_results(i).metric_array(:,1);
    dissipation_vec = master_struct(1).simulation_results(i).energy_per_cycle.*...
        master_struct(1).simulation_results(i).cycle_flux_vec;
    real_flag = master_struct(1).simulation_results(i).real_flag_vec;
    % apply real filter
    non_eq_rates = non_eq_rates(real_flag,:);
    dissipation_vec = dissipation_vec(real_flag);
    information_rate = information_rate(real_flag);
    fidelity_vec = fidelity_vec(real_flag);
    % perform PCA of rates
    [coeff, score_mat] = pca(non_eq_rates);  
    % generate normalize dissipation vec
    % make scatter plot 0f 1st two components
    % get top 100 networks for information and fidelity metrics
    [~, siF] = sort(fidelity_vec,'descend');
    [~, siI] = sort(information_rate,'descend');
    
    pca_scatter = figure;
    hold on
    colormap(cmap);
    scatter(score_mat(:,1),score_mat(:,2),10,[.5 .5 .5],'filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',0)
    scatter(score_mat(siF<=100,1),score_mat(siF<=100,2),20,cmap(1,:),'filled','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',0)
    scatter(score_mat(siI<=100,1),score_mat(siI<=100,2),20,cmap(end,:),'filled','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',0)
end
    
    