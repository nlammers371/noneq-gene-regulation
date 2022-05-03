% script to call core parameter sweep function to examine tradeoffs between
% different network behaviors

clear 
close all
addpath(genpath('../utilities/'))

currentPath = pwd;
if strcmp(currentPath(1:7),'P:\Nick')
    DropboxPath = 'S:\Nick\Dropbox\Nonequilibrium\Nick\SweepOutput\';
else    
    DropboxPath = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\SweepOutput\';
end    
    
FigPath = [DropboxPath 'multi_ns_site_sensitivities' filesep];
mkdir(FigPath);

[~,metric_names] = calculateMetricsNumeric([]);
nStateVec = [6 18 54 162];

rate_index = find(strcmp(metric_names,'Production Rate'));
spec_index = find(strcmp(metric_names,'Specificity'));
spec_alt_index = find(strcmp(metric_names,'specFactorAlt'));
sharp_right_index = find(strcmp(metric_names,'SharpnessRight'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
sharpness_norm_index = find(strcmp(metric_names,'SharpnessRightNorm'));
ir_index = find(strcmp(metric_names,'DecisionRateNorm'));
decision_time_index = find(strcmp(metric_names,'DecisionTimeNorm'));
phi_index = find(strcmp(metric_names,'Phi'));
cw_index = find(strcmp(metric_names,'CW'));

% set path to results
ReadPath = [DropboxPath 'parameter_sweeps_multi_wrong_site/'];

% load results
results_struct = struct;
for n = 1:length(nStateVec)-1
    load([ReadPath 'ir_sharp_spec_n'  num2str(nStateVec(n)) '.mat'])
    fnames = fieldnames(results_struct_ir);
    for f = 1:length(fnames)
        results_struct(n).(fnames{f}) = results_struct_ir.(fnames{f});
    end
end    

%% calculate upper bounds
close all
cw_bins = linspace(0,6,26)';
cw_axis = cw_bins(1:end-1)+diff(cw_bins)/2;
% initialize arrays
ir_global_max_array = NaN(length(cw_axis),length(results_struct));
sharp_ir_max_array = NaN(length(cw_axis),length(results_struct));
spec_ir_max_array = NaN(length(cw_axis),length(results_struct));

ir_sharp_max_array = NaN(length(cw_axis),length(results_struct));
sharp_sharp_max_array = NaN(length(cw_axis),length(results_struct));
ir_spec_max_array = NaN(length(cw_axis),length(results_struct));
spec_spec_max_array = NaN(length(cw_axis),length(results_struct));

for r = 1:length(results_struct)
    % global optimum
    cw_vec_ir = [];
    ir_vec_ir = [];
    sh_vec_ir = [];
    sp_vec_ir = [];
    r_vec_ir = [];
    for i = 1:length(results_struct(r).sim_struct_ir)
        cw_vec_ir = vertcat(cw_vec_ir,results_struct(r).sim_struct_ir(i).metric_array(:,cw_index));
        ir_vec_ir = vertcat(ir_vec_ir,results_struct(r).sim_struct_ir(i).metric_array(:,ir_index));
        sh_vec_ir = vertcat(sh_vec_ir,results_struct(r).sim_struct_ir(i).metric_array(:,sharpness_index));
        sp_vec_ir = vertcat(sp_vec_ir,results_struct(r).sim_struct_ir(i).metric_array(:,spec_index));
        r_vec_ir = vertcat(r_vec_ir,results_struct(r).sim_struct_ir(i).metric_array(:,rate_index));
    end
    
    cw_groups_ir = discretize(cw_vec_ir,cw_bins);
    for c = 1:length(cw_axis)
        [ir_global_max_array(c,r),ir_i] = nanmax(ir_vec_ir.*(1*(cw_groups_ir==c&sh_vec_ir>=0)));
        sharp_ir_max_array(c,r) = sh_vec_ir(ir_i) ./ (r_vec_ir(ir_i)*(1-r_vec_ir(ir_i)));
        spec_ir_max_array(c,r) = sp_vec_ir(ir_i) + cw_vec_ir(ir_i)*(r-1);        
    end
    
    % sharpness optimum   
    cw_vec_sh = [];
    ir_vec_sh = [];
    sh_vec_sh = [];
    sp_vec_sh = [];
    r_vec_sh = [];
    for i = 1:length(results_struct(r).sim_struct_ir)
        cw_vec_sh = vertcat(cw_vec_sh,results_struct(r).sim_struct_sh(i).metric_array(:,cw_index));
        ir_vec_sh = vertcat(ir_vec_sh,results_struct(r).sim_struct_sh(i).metric_array(:,ir_index));
        sh_vec_sh = vertcat(sh_vec_sh,results_struct(r).sim_struct_sh(i).metric_array(:,sharpness_index));
        sp_vec_sh = vertcat(sp_vec_sh,results_struct(r).sim_struct_sh(i).metric_array(:,spec_index) + ...
                    results_struct(r).sim_struct_sh(i).metric_array(:,cw_index).^(r-1));
                  
        r_vec_sh = vertcat(r_vec_sh,results_struct(r).sim_struct_sh(i).metric_array(:,rate_index));
    end
% 
    sh_vec_sh_norm = sh_vec_sh.*r_vec_sh.*(1-r_vec_sh);
    % filter for only those networks which are sufficiently sharp    
    cw_groups_sh = discretize(cw_vec_sh,cw_bins);
    for c = 1:length(cw_axis)
        group_filter = cw_groups_sh==c & sh_vec_sh>=0 & ~isnan(ir_vec_sh);
        if any(group_filter)
            sh_filt = sh_vec_sh_norm(group_filter);
            ir_filt = ir_vec_sh(group_filter);
            sh_val = prctile(sh_filt,95);
            ir_sharp_max_array(c,r) = nanmax(ir_filt(sh_filt>=sh_val));            
        end
    end
%     
%     % specificity optimum   
    cw_vec_sp = [];
    ir_vec_sp = [];
    sh_vec_sp = [];
    sp_vec_sp = [];
    for i = 1:length(results_struct(r).sim_struct_ir)
        cw_vec_sp = vertcat(cw_vec_sp,results_struct(r).sim_struct_sp(i).metric_array(:,cw_index));
        ir_vec_sp = vertcat(ir_vec_sp,results_struct(r).sim_struct_sp(i).metric_array(:,ir_index));
        sh_vec_sp = vertcat(sh_vec_sp,results_struct(r).sim_struct_sp(i).metric_array(:,sharpness_norm_index));
        sp_vec_sp = vertcat(sp_vec_sp,results_struct(r).sim_struct_sp(i).metric_array(:,spec_index) + ...
                    results_struct(r).sim_struct_sp(i).metric_array(:,cw_index).^(r-1));
    end

    % filter for only those networks which are sufficiently sharp    
%     cw_groups_sp = discretize(cw_vec_sp,cw_bins);
%     for c = 1:length(cw_axis)
%         group_filter = cw_groups_sp==c & sh_vec_sp<=1;
%         if any(group_filter)
%             ir_spec_max_array(c,r) = nanmax(ir_vec_sp(group_filter));
%         end
%     end
end    
%%
cv_factor = log2(exp(1));
linetypes = {'-.','--','-'};
global_max_fig = figure;
cmap = brewermap([],'Set2');
hold on
p = [];
for r = length(results_struct):-1:1
    p(r) = plot(10.^cw_axis, ir_global_max_array(:,r)*cv_factor,linetypes{r},'Color',cmap(r,:),'LineWidth',2);
end
set(gca,'yscale','log');
set(gca,'xscale','log');
grid on
ylim([1e-6 1e-1])
set(gca,'Fontsize',14)

xlabel('c_w/c_r')
ylabel('IR (bits per cycle)')
legend(p,'1 spec. site','1 spec. 1 non-spec.','1 spec. 2 non-spec.')

set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
set(gca,'xtick',[1 1e2 1e4 1e6],'xticklabels',{'\alpha^0','\alpha^1','\alpha^2','\alpha^3'})

global_max_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
saveas(global_max_fig,[FigPath 'global_ir_vs_cw.png'])