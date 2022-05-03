clear 
close all
addpath(genpath('../utilities/'))

% set basic parameters
nStates = 6;
rate_bounds = repmat([-6 ; 6],1,3*nStates-4); % constrain transition rate magnitude
[~,metric_names] = calculateMetricsMultiState([]);

% make sure we're linked to the appropriate function subfolder
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(['../utilities/metricFunctions/n' num2str(nStates) '_OR/']));
% addpath(genpath('P:\Nick\projects\noneq-transcription\src\utilities\metricFunctions\n6_OR\from4StateScripts\'));

% define save path
PWD = pwd;
if strcmpi(PWD(1),'P')
    DropboxFolder = 'S:\Nick\Dropbox\Nonequilibrium\Nick\manuscript\';
else
    DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\manuscript\';
end

FigPath = [DropboxFolder 'info_titration' filesep];
mkdir(FigPath);         

% get index of useful metrics
flux_index = find(strcmp(metric_names,'Flux'));
rate_index = find(strcmp(metric_names,'Production Rate'));
spec_index = find(strcmp(metric_names,'Specificity'));
spec_alt_index = find(strcmp(metric_names,'specFactorAlt'));
sharp_right_index = find(strcmp(metric_names,'SharpnessRight'));
sharp_right_norm_index = find(strcmp(metric_names,'SharpnessRightNorm'));
decision_rate_index = find(strcmp(metric_names,'DecisionRateNorm'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
precision_index = find(strcmp(metric_names,'Precision'));
precision_right_index = find(strcmp(metric_names,'PrecisionRight'));

% set sim options
sweep_options = {'n_seeds',5,'n_iters_max',50,'nStates',nStates};

% calculate sensitivity bound
beta = 100;
f0_vec = logspace(log10(beta),log10(beta^2));
seq = 1/4;

% specify plot options 
n_plot = 3e3; % number of points to plot
markerAlpha = 0.5; % marker transparency
markerSize = 75; % marker size

%% %% Track info rate vs mu for f0 and s0-optimized networks %%%%

nBins = 101;
if ~exist([FigPath 'info_s0_f0_motif_struct.mat'])    
    info_s0_f0_motif_struct = generate_info_vs_s0_vs_f0_vs_mu_scatters(nBins,sweep_options,beta,FigPath);
else
    load([FigPath 'info_s0_f0_motif_struct.mat'],'info_s0_f0_motif_struct')
end

%%%%%%%%%%%%%%%%%
% generate predicted curves for the two motifs
c_factor = (0.1/1.05)^2;
mu_vec = [info_s0_f0_motif_struct.mu_vec];
x_axis = mu_vec/beta;
x_lim = ([1e-2 1e5]);
x_tick = [1e-2 1e0 1e2 1e4 1e6];
%%%%%%%%%%%%%%%%%%%%%%%%
% make figure
close all
motif_fig = figure;
hold on
cmap = brewermap([],'Set2');

% sweep results for sharpness-optimized
% s1 = scatter(x_axis, [info_s0_f0_motif_struct.info_vec_spec]/c_factor,...
%           0.5*markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(4,:));
% prediction for sharpness-optimized
p3 = plot(x_axis,[info_s0_f0_motif_struct.info_max_vec]/c_factor,'-','Color',[cmap(2,:) 0.5],'LineWidth',2);

i_spec_vec = [info_s0_f0_motif_struct.info_vec_spec]/c_factor;
p1 = plot(x_axis,imgaussfilt(i_spec_vec,1),'--','Color',brighten(cmap(4,:),-0.5),'LineWidth',2);
% 
% sweep for f0-optimized
x_ft = true(size(x_axis));%round(x_axis,2)~=5.37;% & x_axis<=10^2;
i_sharp_vec = [info_s0_f0_motif_struct.info_vec_sharp]/c_factor;

% s2 = scatter(x_axis(x_ft), i_s_vec(x_ft)/c_factor,...
%           0.5*markerSize,'s','MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(5,:));        
               
% f0 prediction
p2 = plot(x_axis(x_ft),imgaussfilt(i_sharp_vec(x_ft),2),'--','Color',brighten(cmap(5,:),-0.5),'LineWidth',2);        


% s3 = scatter(x_axis, [info_s0_f0_motif_struct.info_max_vec]/c_factor,...
%           0.5*markerSize,'d','MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(9,:));            

set(gca,'xscale','log')
set(gca,'yscale','log')
ylim([10^-10 10^0.5])
xlim(x_lim)
xlabel('wrong factor ratio (\mu/\beta)');
ylabel('information rate (V*)')

legend([p2 p1 p3],'s_0 optimized','f_0 optimized','global optimum (neq)','Location','southwest')

grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
set(gca,'xtick',x_tick)
%set(gca,'xticklabels',{'\beta^{0}','\beta^{1}','\beta^{2}','\beta^{3}'})
motif_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(motif_fig,[FigPath 'motif_information_plot.png'])
saveas(motif_fig,[FigPath 'motif_information_plot.pdf'])

 %%%%%% plot max eq and neq info rate 
close all
info_mu_fig = figure;
hold on
                      
% scatter(x_axis, [info_s0_f0_motif_struct.info_max_vec]/c_factor,...
%           0.5*markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(2,:)); 
plot(x_axis, [info_s0_f0_motif_struct.info_max_vec]/c_factor,...
           'Color',cmap(2,:),'LineWidth',2);
% scatter(x_axis, [info_s0_f0_motif_struct.info_max_vec_eq]/c_factor,...
%           0.5*markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(3,:));         
plot(x_axis, [info_s0_f0_motif_struct.info_max_vec_eq]/c_factor,...
           'Color',cmap(3,:),'LineWidth',2);

set(gca,'xscale','log')
set(gca,'yscale','log')
ylim([10^-10 10^0.5])
xlim(x_lim)
xlabel('wrong factor ratio (\mu/\beta)');
ylabel('information rate (V*)')

% yyaxis right
% plot(x_axis, [info_s0_f0_motif_struct.info_max_vec]./[info_s0_f0_motif_struct.info_max_vec_eq],'--',...
%            'Color',cmap(8,:),'LineWidth',2);

set(gca,'yscale','log')

legend('non-equilibrium','equilibrium','location','southwest')

grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
set(gca,'xtick',x_tick)
%set(gca,'xticklabels',{'\beta^{0}','\beta^{1}','\beta^{2}','\beta^{3}'})
info_mu_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(info_mu_fig,[FigPath 'information_vs_mu_plot.png'])
saveas(info_mu_fig,[FigPath 'information_vs_mu_plot.pdf'])
%%
f_factor_vec = 10.^([info_s0_f0_motif_struct.f0_vec_info])*beta;
seq2 = 0.25;

f01 = beta;
f02 = beta^2;
s_bound_1 = (f01./mu_vec) ./ (1 + f01./mu_vec) .* seq2 ; 
s_bound_2 = (f_factor_vec./mu_vec) ./ (1 + f_factor_vec./mu_vec) .* seq2 .* (beta-1)./beta .* (beta + f_factor_vec) ./ (f_factor_vec- 1); 
sharp_pd_vec = s_bound_2./s_bound_1;

info_mu_ratio_fig = figure;
hold on
                      
% scatter(x_axis, [info_s0_f0_motif_struct.info_max_vec]./[info_s0_f0_motif_struct.info_max_vec_eq],...
%           0.5*markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(8,:));     
plot(x_axis, imgaussfilt([info_s0_f0_motif_struct.info_max_vec]./[info_s0_f0_motif_struct.info_max_vec_eq],2),...
           'Color','k','LineWidth',2);
         
plot(x_axis, imgaussfilt(sharp_pd_vec.^2,2), '--','Color','k','LineWidth',1.5);         

set(gca,'xscale','log')
set(gca,'yscale','log')

xlabel('wrong factor ratio (\mu/\beta)');
ylabel('information gain (g)')

legend('numerical results','sharpness-based prediction','Location','southeast')

grid on
set(gca,'FontSize',14)
set(gca,'xtick',x_tick)
%set(gca,'xticklabels',{'\beta^{0}','\beta^{1}','\beta^{2}','\beta^{3}'})
set(gca,'Color',[228,221,209]/255) 
xlim(x_lim)
ylim([10^0 10^4.5])
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

info_mu_ratio_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(info_mu_ratio_fig,[FigPath 'noneq_gain_vs_mu_plot.png'])
saveas(info_mu_ratio_fig,[FigPath 'noneq_gain_vs_mu_plot.pdf'])

%%

% define performance limits
prec_eq = 8;
prec_neq = 16;

s0_eq = 1;
s0_neq = 4;

f0_eq = 1;
f0_neq = beta;

% extract vectors
s_factor_vec = ([info_s0_f0_motif_struct.s0_vec_info]).^2;
f_factor_vec = 10.^([info_s0_f0_motif_struct.f0_vec_info]);
p_factor_vec = exp([info_s0_f0_motif_struct.prec_vec_info]).^2;

close all

metric_contribution_fig = figure;
hold on
         
s1 = scatter(x_axis, imgaussfilt((s_factor_vec-s0_eq)/(s0_neq-s0_eq),2),...
  0.65*markerSize,'s','MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',1, 'MarkerFaceColor',cmap(5,:)); 

s2 = scatter(x_axis, imgaussfilt((f_factor_vec-f0_eq)/(f0_neq-f0_eq),2),...
  0.65*markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',1, 'MarkerFaceColor',cmap(4,:));     
% 
s3 = scatter(x_axis, imgaussfilt((p_factor_vec-prec_eq)/(prec_neq-prec_eq),2),...
  0.65*markerSize,'^','MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',1, 'MarkerFaceColor',cmap(6,:));     



set(gca,'xscale','log')

xlabel('wrong factor ratio (\mu/\beta)');
ylabel('optimization score')

legend([s1 s2 s3],'top. sharpness (s_0)','activator fidelity (f_0)','precision (1/\sigma_r)','Location','southwest')

grid on
set(gca,'FontSize',14)
% set(gca,'xtick',[1e-4 1e-2 1e0 1e2 1e4])
set(gca,'xtick',x_tick)
%set(gca,'xticklabels',{'\beta^{0}','\beta^{1}','\beta^{2}','\beta^{3}'})
set(gca,'Color',[228,221,209]/255) 
xlim(x_lim)
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

metric_contribution_fig.InvertHardcopy = 'off';
set(gcf,'color','w');


saveas(metric_contribution_fig,[FigPath 'noneq_metric_contributions_vs_mu_plot.png'])
saveas(metric_contribution_fig,[FigPath 'noneq_metric_contributions_vs_mu_plot.pdf'])
% hold on
% scatter(x_axis, )

% Look at predicted decision times
% set parameters
c1 = 1.1;
c0 = 1;
epsilon = 0.32;


info_rates_eq = vertcat(info_s0_f0_motif_struct.info_max_rates_eq);
info_rates_neq = vertcat(info_s0_f0_motif_struct.info_max_rates);

% generate input cell arrays
valMatEq = [c0*ones(size(info_rates_eq,1),1) info_rates_eq];
valCellEq = mat2cell(valMatEq,size(valMatEq,1),ones(1,size(valMatEq,2)));

valMatNeq = [c0*ones(size(info_rates_neq,1),1) info_rates_neq];
valCellNeq = mat2cell(valMatNeq,size(valMatNeq,1),ones(1,size(valMatNeq,2)));

% calculate cycle times
TauOnEq = TauONFunction(valCellEq{:});
TauOffEq = TauOFFFunction(valCellEq{:});
TauCycleEq = TauOffEq+TauOnEq;

TauOnNeq = TauONFunction(valCellNeq{:});
TauOffNeq = TauOFFFunction(valCellNeq{:});
TauCycleNeq = TauOffNeq+TauOnNeq;

% normalize to standardize cycle times
info_rates_eq_norm = info_rates_eq.*TauCycleEq/60;
info_rates_neq_norm = info_rates_neq.*TauCycleNeq/60;

% calculate decision times for a range of error tolerance
[~,tau_vec_eq,~] = calculateDecisionMetrics(info_rates_eq_norm,c1,c0,epsilon);
[~,tau_vec_neq,~] = calculateDecisionMetrics(info_rates_neq_norm,c1,c0,epsilon);

close all

decision_time_fig = figure;

hold on
plot(x_axis,tau_vec_eq/60,'Color',cmap(3,:),'LineWidth',3)
plot(x_axis,tau_vec_neq/60,'Color',cmap(2,:),'LineWidth',3)

legend('non-equilibrium','equilibrium','Location','northwest')
set(gca,'FontSize',14)
% set(gca,'xtick',[1e-4 1e-2 1e0 1e2 1e4])
set(gca,'Color',[228,221,209]/255) 

xlabel('wrong factor ratio (\mu/\beta)');
ylabel('decision time (minutes)')

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

decision_time_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
xlim(xlim)
set(gca,'xtick',x_tick)
%set(gca,'xticklabels',{'\beta^{0}','\beta^{1}','\beta^{2}','\beta^{3}'})
set(gca,'ytick',[1e2 1e4 1e6 1e8 1e10])
ylim([1e1 10^10])
xlim(x_lim)
set(gca,'yscale','log')
set(gca,'xscale','log')

grid on

saveas(decision_time_fig,[FigPath 'decision_time_vs_mu.png'])
saveas(decision_time_fig,[FigPath 'decision_time_vs_mu.pdf'])

%%
% scatter(x_axis, [info_s0_f0_motif_struct.i_vec_spec]/c_factor,...
%           0.5*markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(4,:));

r_vec = [info_s0_f0_motif_struct.rate_i_max_vec];
close all
figure;
hold on
scatter(x_axis,2*([info_s0_f0_motif_struct.i_max_vec]/c_factor) ./ [info_s0_f0_motif_struct.sharp_i_max_vec].^2 ,'filled'); 
% scatter(x_axis,[info_s0_f0_motif_struct.spec_i_max_vec],'filled'); 
% scatter(x_axis,s_bound_fun(beta^2).^2*16)
set(gca,'xscale','log')
set(gca, 'XDir','reverse')
%%
[sim_info_neq, sim_struct_neq] = param_sweep_multi([spec_index precision_index],sweep_options{:},...
                                                  'half_max_flag',true,'wrongFactorConcentration',1e2,'equilibrium_flag',...
                                                  false,'n_sim',5,'useParpool',0,'specFactor',beta);

%%                                                
results_array = vertcat(sim_struct_neq.metric_array);
% close all;
% figure;
hold on
scatter(results_array(:,spec_index),exp(results_array(:,precision_index)).^2)  
% ylim([10 20])