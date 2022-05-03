clear 
close all
addpath('../utilities')
dataPath = '../../out/param_optimization/';
figPath = '../../fig/param_optimization/';
mkdir(figPath);
% load data
activator_flag = 1;
load([dataPath 'tau_fold_dependence_act' num2str(activator_flag) '.mat'],'sim_struct');
% plot decision times vs. fold difference
c_profile = sim_struct.c_profile;
x_axis = c_profile(1)./c_profile(2:end);
decision_time_array = sim_struct.decision_time_array;
n_models = size(decision_time_array,3);

dt_fig = figure;
cmap = flipud(brewermap(9,'Set2'));
color_indices = fliplr([2,7,8,5]);
hold on
for n = 1:n_models
    if n == 1
        lt = '--';
        cp = 'black';
    else
        lt = '-';
        cp = cmap(color_indices(n-1),:);
    end
    slice = decision_time_array(:,:,n);
    dt_vec = nanmin(slice,[],2);
    plot(x_axis,dt_vec,lt,'Color',cp,'LineWidth',1.5);
end
legend('equilibrium','neq (specific)', 'neq (non-specific)', 'neq (full)')
xlabel('fold concentration difference')
ylabel('decision time (seconds)')
set(gca,'Fontsize',14,'YScale','log','XScale','log')
grid on
saveas(dt_fig,[figPath 'fold_diff_d_times_act' num2str(activator_flag) '.png'])