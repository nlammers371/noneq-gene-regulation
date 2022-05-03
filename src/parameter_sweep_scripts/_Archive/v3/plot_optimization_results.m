% Script to plot results of point optimization routines for parameters
clear
close all

readPath = '../../out/param_sweep/';
figPath = '../../fig/param_sweep/';
mkdir(figPath);

load([readPath 'param_point_optimization.mat'])
metric_names = {'fidelity','dynamicRange','info_rate'};
cmap = brewermap(128,'Spectral');
for i = 1:numel(opt_results)
    activator_flag = i == 1;
    metric_array = opt_results(i).metric_array;
    fid_results = reshape(nanmax(metric_array(:,1,:)),1,[]);
    range_results = reshape(nanmax(metric_array(:,2,:)),1,[]);
    info_results = reshape(nanmax(metric_array(:,3,:)),1,[]);
    index_vec = 1:numel(fid_results);
    
    % fidelity vs. sharpness
    f1 = figure;
    colormap(cmap)
    s = scatter(fid_results+rand(size(fid_results))*.005,range_results+rand(size(fid_results))*.005...
        ,40,index_vec,'filled','MarkerEdgeColor','black');
    s.MarkerFaceAlpha = .4;
    h = colorbar;
    ylabel('sharpness')
    xlabel('fidelity')
    ylabel(h,'model ID')
    set(gca,'Fontsize',14)
    grid on
    ylim([0,.55])
    xlim([0,2])
    saveas(f1,[figPath 'fid_v_sharp_opt_results_act' num2str(activator_flag) '.tif'])
    
    % sharpness vs. info
    f2 = figure;
    colormap(cmap)
    s = scatter(fid_results+rand(size(fid_results))*.005,info_results+rand(size(fid_results))*.005...
        ,40,index_vec,'filled','MarkerEdgeColor','black');
    s.MarkerFaceAlpha = .4;
    h = colorbar;
    xlabel('sharpness')
    ylabel('information rate')
    ylabel(h,'model ID')
    set(gca,'Fontsize',14)
    grid on
    ylim([0,6])
    saveas(f2,[figPath 'sharp_v_inf_opt_results_act' num2str(activator_flag) '.tif'])
end
