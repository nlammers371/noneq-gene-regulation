clear
close all
% specify paths
addpath('../utilities')
readPath = '../../out/parameter_search_v2/';
figPath = '../../fig/parameter_search_v2/';
% suffix_cell = {'dmdc_InvSigma', 'pointwise-error_average-error', 'fidelity_dynamicRange',...
%     'fidelity_info-rate'};
% x_names = {'sharpness',  'pointwise precision', 'fidelity','fidelity'};
% y_names = {'\sigma_m^{-1}',  'average precision', 'dynamic range', 'information rate'};
suffix_cell = {'fold-concentration_fold-production','dmdc_InvSigma','decisionRate_boundPosition'};
x_names = {'log(c_1/c_0)', 'sharpness', 'decision rate'};
y_names = {'log(m_1/m_0)', 'precision (\sigma_m^{-1})', 'boundary position',};

mkdir(figPath)
metric_names = {'pointwise-error','average-error','dynamicRange','fidelity','info-rate','dmdc',...
    'InvSigma','decisionRate','boundPosition','fold-concentration','fold-production'};
% load data
master_struct = struct;
for i = 1:numel(suffix_cell)
    load([readPath 'edge_sim_results_' suffix_cell{i} '.mat'])
    fnames = fieldnames(simulation_results);
    master_struct(i).simulation_results = simulation_results;
    master_struct(i).type = suffix_cell{i};
end


type_names = {'activator', 'repressor'};
% plot metric dispersion
cmap1 = brewermap(9,'Set2');
cmap2 = brewermap(128,'Spectral');
mSize = 15;
for j = 3%1:numel(master_struct)
    name = suffix_cell{j};
    name1 = name(1:strfind(name,'_')-1);
    name2 = name(strfind(name,'_')+1:end);
    ind1 = find(ismember(metric_names, name1));
    ind2 = find(ismember(metric_names, name2));
    for i = 1:numel(simulation_results)    
        if i == 1
            cn = cmap2(115,:);
            ce = cmap2(90,:);
        else
            cn = cmap1(3,:);
            ce = cmap1(2,:);
        end
        % extract relevant results
        neq1 = reshape(master_struct(j).simulation_results(i).metric_array_noneq(:,ind1,:,:),[],1);
        neq2 = reshape(master_struct(j).simulation_results(i).metric_array_noneq(:,ind2,:,:),[],1);
        eq1 = reshape(master_struct(j).simulation_results(i).metric_array_eq(:,ind1,:,:),[],1);
        eq2 = reshape(master_struct(j).simulation_results(i).metric_array_eq(:,ind2,:,:),[],1);        
        if ind2 == 7
            neq2 = exp(neq2);
            eq2 = exp(eq2);
        end
        % (1) plot inverse noise vs sharpness eq x neq
        noise_vs_sharpness = figure;
        hold on           
        scatter(neq1,neq2,mSize,'MarkerFaceColor',cn,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.5);
        scatter(eq1,eq2,mSize,'MarkerFaceColor',ce,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.5);
        xlabel(x_names{j})
        ylabel(y_names{j})
        legend('non-equilibrium','equilibrium','Location','best')
        grid on
        set(gca,'Fontsize',12)
        if ind1 == 6
            xlim([.05 .5])   
        end
        saveas(noise_vs_sharpness,[figPath type_names{i} '_' name '.tif'])
    end
end
    
%     
%         % (2) plot inverse noise vs sharpness eq 
%         eq_noise_vs_sharpness = figure;
%         colormap(flipud(cmap))
%         hold on           
%         scatter(eq_sharpness,eq_inv_sigma,mSize,eq_sharpness.*eq_inv_sigma,'filled','MarkerEdgeAlpha',0,'MarkerFaceAlpha',.5);
%         xlabel('sharpness')
%         ylabel('\sigma_m^{-1}')
%         h = colorbar;
%         ylabel(h,'\sigma_c^{-1}')    
%         grid on
%         set(gca,'Fontsize',12)
%         xlim([.05 .45])    
%     %     ylim([0 1])
%         saveas(eq_noise_vs_sharpness,[figPath model_names{i} '_noise_v_sharpness_eq.tif'])
%     %     saveas(eq_noise_vs_sharpness,[figPath model_names{i} '_noise_v_sharpness_eq.pdf'])
%     %     
%         % (3) plot inverse noise vs sharpness neq 
%         neq_noise_vs_sharpness = figure;
%         colormap(flipud(cmap))
%         hold on           
%         scatter(neq_sharpness,neq_inv_sigma,mSize,neq_sharpness.*neq_inv_sigma,'filled','MarkerEdgeAlpha',0,'MarkerFaceAlpha',.5);
%         xlabel('sharpness')
%         ylabel('\sigma_m^{-1}')
%         h = colorbar;
%         ylabel(h,'\sigma_c^{-1}')    
%         grid on
%         set(gca,'Fontsize',12)
%         xlim([.05 .45])    
%     %     ylim([0 1])
%         saveas(neq_noise_vs_sharpness,[figPath model_names{i} '_noise_v_sharpness_neq.tif'])    
% 
%         % extract vectors
%         neq_fidelity = reshape(master_struct(2).simulation_results(i).metric_array_noneq(:,7,:),[],1);
%         neq_info = reshape(master_struct(2).simulation_results(i).metric_array_noneq(:,8,:),[],1);
%         eq_fidelity = reshape(master_struct(2).simulation_results(i).metric_array_eq(:,7,:),[],1);
%         eq_info = reshape(master_struct(2).simulation_results(i).metric_array_eq(:,8,:),[],1);
% 
%         % fidelity vs decision rate
%         fid_v_info_fig = figure;
%         hold on
%         scatter(neq_fidelity,neq_info,mSize,'MarkerFaceColor',cmap(110,:),'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.5);
%         scatter(eq_fidelity,eq_info,mSize,'MarkerFaceColor',cmap(10,:),'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.5);
%         xlabel('fidelity')
%         ylabel('information rate')
%         grid on
%         set(gca,'Fontsize',12)
%         legend('non-equilibrium','equilibrium','Location','northwest')
%     %     axis([0 2 0 15])
%         saveas(fid_v_info_fig,[figPath  model_names{i} '_fidelity_v_info.tif'])
%     %     saveas(fid_v_info_fig,[figPath  model_names{i} '_fidelity_v_info.pdf'])
% 
% 
%     end
% end   