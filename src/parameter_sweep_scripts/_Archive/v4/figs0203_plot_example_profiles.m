% script to make plots of example profiles
clear
close all
addpath('../utilities')
readPath = '../../out/param_sweep/';
figPath = '../../fig/example_profiles/';
mkdir(figPath);
suffix_cell = {'fidelity_dynamicRange'};

% first look at high-performing fidelity network
c_range = linspace(2/11,20/11,101);
suffix_ind = 1;
mci = 1;

load([readPath 'param_sweep_results_' suffix_cell{suffix_ind} '.mat'])
%%% extract rate and metric arrays
eq_ind = 1;
neq_ind = 6;
metric_indices = [3 4];
rate_array_neq = reshape(permute(simulation_results(2).rate_array(:,:,:,:,neq_ind),[1 4 3 2]),[],8);
metric_array_neq = reshape(permute(simulation_results(2).metric_array(:,metric_indices,:,:,neq_ind),[1 4 3 2]),[],2);           
rate_array_eq = reshape(permute(simulation_results(2).rate_array(:,:,:,:,eq_ind),[1 4 3 2]),[],8);
metric_array_eq = reshape(permute(simulation_results(2).metric_array(:,metric_indices,:,:,eq_ind),[1 4 3 2]),[],2);           
% find high-performing networl
[max_f, max_i_neq] = max(metric_array_neq(:,1));
[~, max_i_eq] = max(metric_array_eq(:,1));
rates_neq = rate_array_neq(max_i_neq,:);
rates_eq = rate_array_eq(max_i_eq,:);
% generate profiles for eq and noneq 
[~, profile_neq, dd_prof_neq, var_vec_neq] = calculateMetricsConstrained(rates_neq,c_range,0);        
[~, profile_eq, dd_prof_eq, var_vec_eq] = calculateMetricsConstrained(rates_eq,c_range,0);        
% plot
p_fold_neq = profile_neq ./ profile_neq(end);
p_fold_eq = profile_eq ./ profile_eq(end);

fidelity_fold_fig = figure;
hold on
plot(c_range, p_fold_neq);
plot(c_range, p_fold_eq);
scatter(c_range([1 end 1 end]),[p_fold_neq([1 end]) p_fold_eq([1 end])],'MarkerFaceColor','black',...
    'MarkerEdgeAlpha',0)
grid on
xlabel('[R]')
ylabel('normalized transcription rate')
legend('non-eq','eq')
set(gca,'Fontsize',14)
saveas(fidelity_fold_fig,[figPath 'fidelity_fold_profile.tif'])

fidelity_abs_fig = figure;
hold on
plot(c_range, profile_neq);
plot(c_range, profile_eq);
scatter(c_range([1 end 1 end]),[profile_neq([1 end]) profile_eq([1 end])],'MarkerFaceColor','black',...
    'MarkerEdgeAlpha',0)
grid on
xlabel('[R]')
ylabel('transcription rate')
legend('non-eq','eq')
set(gca,'Fontsize',14)
% ylim([0 1])
saveas(fidelity_abs_fig,[figPath 'fidelity_abs_profile.tif'])
%%%
% Now dynamic range
[max_d_neq, max_i_neq] = max(metric_array_neq(:,2));
[max_d_eq, max_i_eq] = max(metric_array_eq(:,2));
dr_rates_neq = rate_array_neq(max_i_neq,:);
dr_rates_eq = rate_array_eq(max_i_eq,:);

[~, dr_profile_neq, dd_prof, var_vec] = calculateMetricsConstrained(dr_rates_neq,c_range,0);        
[~, dr_profile_eq, ~, ~] = calculateMetricsConstrained(dr_rates_eq,c_range,0);        

drange_fig = figure;
hold on
plot(c_range, dr_profile_neq);
plot(c_range, dr_profile_eq);
scatter(c_range([1 end 1 end]),[dr_profile_neq([1 end]) dr_profile_eq([1 end])],'MarkerFaceColor','black',...
    'MarkerEdgeAlpha',0)
grid on
xlabel('[R]')
ylabel('transcription rate')
legend('non-eq','eq')
set(gca,'Fontsize',14)
saveas(drange_fig,[figPath 'dynamic_range_profile.tif'])

%%% Look at network performance as amplifier

[~, dr_profile_amp_neq, dd_prof, var_vec] = calculateMetricsConstrained(dr_rates_neq,c_range.^2,0);        
[~, dr_profile_amp_eq, ~, ~] = calculateMetricsConstrained(dr_rates_eq,c_range.^2,0);        
drange_fig = figure;
hold on
% plot(c_range, dr_profile_neq,'Color',cmap(1,:));
plot(c_range, dr_profile_amp_neq,'LineWidth',1.5);
% plot(c_range, dr_profile_eq,'Color',cmap(3,:));
plot(c_range, dr_profile_amp_eq,'LineWidth',1.5);

grid on
xlabel('[R]')
ylabel('transcription rate')
legend('non-eq','eq')
set(gca,'Fontsize',14)
saveas(drange_fig,[figPath 'dynamic_range_amp_profile1.png'])

drange_fig = figure;
hold on
% plot(c_range, dr_profile_neq,'Color',cmap(1,:));
plot(c_range, dr_profile_neq,'LineWidth',1.5);
% plot(c_range, dr_profile_eq,'Color',cmap(3,:));
plot(c_range, dr_profile_eq,'LineWidth',1.5);

grid on
xlabel('[R]')
ylabel('transcription rate')
legend('non-eq','eq')
set(gca,'Fontsize',14)
saveas(drange_fig,[figPath 'dynamic_range_amp_profile2.png'])

%% illustrative noise plot
cmap = brewermap(9,'Set2');
ub_p = dr_profile_neq(1:end-1) + sqrt(var_vec);
lb_p = dr_profile_neq(1:end-1) - sqrt(var_vec);

noise_fig = figure;
hold on
f = fill([c_range(1:end-1) fliplr(c_range(1:end-1))],[ub_p fliplr(lb_p)],cmap(2,:));
f.FaceAlpha = .5;
f.EdgeAlpha = 0;
plot(c_range(1:end-1),dr_profile_neq(1:end-1),'Color',cmap(2,:),'LineWidth',2)
grid on
xlim([c_range(1) c_range(end)])
xlabel('[R]')
ylabel('transcription rate')
set(gca,'Fontsize',14)
saveas(noise_fig,[figPath 'noise_profile.tif'])