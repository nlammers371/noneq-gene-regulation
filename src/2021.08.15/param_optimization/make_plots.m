clear
close all
addpath('../utilities')
dataPath = '../../out/param_optimization/';
figPath = '../../fig/param_optimization/';
mkdir(figPath);

% specify burn-in period
f_diff = 30;
% load point-wise data
point_struct = struct;
load([dataPath 'rate_optimization_results_f' num2str(f_diff) '_act1.mat'],'sim_struct');
point_struct(1).sim_struct = sim_struct;
point_struct(1).activator_flag = 1;
load([dataPath 'rate_optimization_results_f' num2str(f_diff) '_act0.mat'],'sim_struct');
point_struct(2).sim_struct = sim_struct;
point_struct(2).activator_flag = 0;
% define colormap
cmap = flipud(brewermap(128,'RdYlBu'));
inc = 20;
for i = 1:numel(point_struct)
    activator_flag = i==1;
    decision_time_array = point_struct(i).sim_struct.decision_time_array(:,1:end,:);
    x_vec = linspace(0,1,size(decision_time_array,1));
    x_long = repmat(x_vec,1,size(decision_time_array,2));
    scatter_fig = figure;
    hold on
    eq_slice = decision_time_array(:,:,1);
    min_eq = nanmin(eq_slice,[],2);        
    plot(x_vec,min_eq/60,'Color','black')
    for j = 2:size(decision_time_array,3)
        x_rand = x_vec + rand(size(x_vec))*.005;
        d_slice = decision_time_array(:,:,j);
        min_slice = nanmin(d_slice,[],2);        
        scatter(x_rand,min_slice/60,20,'MarkerFaceColor',cmap(10+inc*(j-1),:),'MarkerFaceAlpha',.5,'MarkerEdgeColor','black');
    end
    grid on    
    xlabel('position')
    ylabel('decision time (min)')
    set(gca,'Fontsize',14,'yscale','log'z)
    saveas(scatter_fig,[figPath 'decision_time_f' num2str(f_diff) '_act' num2str(activator_flag) '.tif'])
end

%% fold dependence
fold_struct = struct;
load([dataPath 'tau_fold_dependence_act1.mat'],'sim_struct');
fold_struct(1).sim_struct = sim_struct;
fold_struct(1).activator_flag = 1;
load([dataPath 'tau_fold_dependence_act0.mat'],'sim_struct');
fold_struct(2).sim_struct = sim_struct;
fold_struct(2).activator_flag = 0;

cmap = flipud(brewermap(9,'Set2'));
for i = 1:numel(fold_struct)
    activator_flag = i==1;
    decision_time_array = fold_struct(i).sim_struct.decision_time_array;
    c_vec = fold_struct(i).sim_struct.c_profile(1) ./ fold_struct(i).sim_struct.c_profile;
    fold_fig = figure;
    hold on
    eq_slice = decision_time_array(:,:,1);
    
    for j = 2:size(decision_time_array,3)        
        d_slice = decision_time_array(:,:,j);
        min_slice = nanmin(d_slice,[],2);        
        plot(c_vec(2:end)',min_slice,'Color',cmap(j,:),'LineWidth',1.5);
    end
    min_eq = nanmin(eq_slice,[],2);        
    plot(c_vec(2:end),min_eq','--','Color','black','LineWidth',1.5)
    grid on    
    xlabel('fold difference')
    ylabel('decision time (sec)')
    set(gca,'Fontsize',12,'xscale','log','yscale','log')
    saveas(fold_fig,[figPath 'fold_dependence_act' num2str(activator_flag) '.tif'])
end



%%
% load full profile data 
full_prof_struct = struct;
load([dataPath 'rate_optimization_full_profile_act1.mat'],'sim_struct');
full_prof_struct(1).sim_struct = sim_struct;
full_prof_struct(1).activator_flag = 1;
load([dataPath 'rate_optimization_full_profile_act0.mat'],'sim_struct');
full_prof_struct(0).sim_struct = sim_struct;
full_prof_struct(0).activator_flag = 0;