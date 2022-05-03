% script to call core parameter sweep function to examine tradeoffs between
% different network behaviors

clear 
close all
addpath(genpath('../utilities/'))

% set basic parameters
nStates = 6;
rate_bounds = repmat([-6 ; 6],1,3*nStates-4); % constrain transition rate magnitude
[~,metric_names] = calculateMetricsSym([]);

% specify function path
functionPath = ['../utilities/metricFunctions/n' num2str(nStates) '_OR/'];

% make sure we're linked to the appropriate function subfolder% make sure we're linked to the appropriate function subfolder
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));

% define save path

DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\manuscript\';
% DropboxFolder = 'S:\Nick\Dropbox\Nonequilibrium\Nick\manuscript\';

FigPath = [DropboxFolder 'optimality_landscape' filesep];
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
cw_index = find(strcmp(metric_names,'CW'));

% set sim options
sweep_options = {'n_seeds',5,'n_iters_max',50,'nStates',nStates};

% calculate sensitivity bound
alpha_factor = 100;
f0_vec = logspace(log10(alpha_factor),log10(alpha_factor^2));
% seq = 1/4;

% specify plot options 
n_plot = 3e3; % number of points to plot
markerAlpha = 0.5; % marker transparency
markerSize = 75; % marker size

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% s0 vs f0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cw = 1e4;
tic
[sim_info_neq, sim_struct_neq] = param_sweep_multi_v2([sharp_right_norm_index spec_index],functionPath,sweep_options{:},...
                                          'half_max_flag',true,'wrongFactorConcentration',cw,...
                                          'equilibrium_flag',false,'specFactor',alpha_factor);

[sim_info_eq, sim_struct_eq] = param_sweep_multi_v2([sharp_right_norm_index spec_index],functionPath,sweep_options{:},...
                                          'half_max_flag',true,'wrongFactorConcentration',cw,...
                                          'equilibrium_flag',true,'specFactor',alpha_factor);                                        
toc     

close all
seq = 1;
c = 1;
%% cw = 1;
% get predicted tradeoff bound
% s0_bound = seq+(alpha_factor.^2+(-1).*f0_vec).*cw.*((-1).*alpha_factor+f0_vec+((-1)+alpha_factor).*f0_vec.*cw).^( ...
%   -1).*seq;
s0_bound = (c.^(-1)+((-1)+alpha_factor).^(-1).*c.^(-1).*(alpha_factor.^2+(-1).*f0_vec).*f0_vec.^(-1)).*seq;             
% s0_bound = seq.*(1+(alpha_factor.^2+(-1).*f0_vec).*u.*((-1).*alpha_factor+f0_vec+((-1)+alpha_factor).*f0_vec.*u).^(-1));

% generate vectors to plot
results_array = vertcat(sim_struct_neq.metric_array);
f0_scatter_vec = 10.^results_array(:,spec_index);
s0_scatter_vec = results_array(:,sharp_right_norm_index);
plot_options = find(s0_scatter_vec>=0 & f0_scatter_vec >= 1);
plot_indices = randsample(plot_options,min([n_plot,length(plot_options)]),false);

% make figure
s0_f0_fig = figure;
hold on
cmap = brewermap([],'Set2');

scatter(f0_scatter_vec(plot_indices),s0_scatter_vec(plot_indices),...
      markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(2,:)); 
    
plot(f0_vec/alpha_factor,s0_bound,'--','Color','k','LineWidth',3)
% ylim([0 0.5])
xlim([1 alpha_factor])
set(gca,'xscale','log')    
xlabel('fidelity (f^{neq}_0/f^{eq}_0)');
ylabel('topological sharpness (s_0^{neq}/s_0^{eq})')
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

s0_f0_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

% saveas(s0_f0_fig,[FigPath 's0_vs_f0.png'])
% saveas(s0_f0_fig,[FigPath 's0_vs_f0.pdf'])


%% %%%%%%%%%%%%%%%% sharpness vs activator fidelity %%%%%%%%%%%%%%%%%%%%%
close all

% [sim_info_neq, sim_struct_neq] = param_sweep_multi([sharpness_index cw_index],sweep_options{:},...
%                                                   'half_max_flag',true,'wrongFactorConcentration',mu_vec(m),'equilibrium_flag',false);

% calculate sensitivity bound
alpha_factor = sim_info_neq.specFactor;
f0_vec = logspace(log10(alpha_factor),log10(alpha_factor^2));
mu_vec = logspace(log10(1),log10(alpha_factor^2),11);


sharpness_cell = cell(1,length(mu_vec));
s0_cell = cell(1,length(mu_vec));
spec_vec = cell(1,length(mu_vec));

s_f0_mu_titration = struct;

if ~exist([FigPath 's_f0_mu_titration.mat'])
    for  m = 1:length(mu_vec)
        tic
        [sim_info_neq, sim_struct_neq] = param_sweep_multi([sharpness_index spec_index],sweep_options{:},...
                                                  'half_max_flag',true,'wrongFactorConcentration',mu_vec(m),'equilibrium_flag',false);

        result_array = vertcat(sim_struct_neq.metric_array);
        s_f0_mu_titration(m).sharpness = result_array(:,sharpness_index);
        s_f0_mu_titration(m).s0 = result_array(:,sharp_right_norm_index);
        s_f0_mu_titration(m).f0 = result_array(:,spec_index);
        toc 
    end
    
    save([FigPath 's_f0_mu_titration.mat'],'s_f0_mu_titration')
else
    load([FigPath 's_f0_mu_titration.mat'],'s_f0_mu_titration')
end  

%% make figure
seq2 = 1/4;
close all 
c = 1;

% make figure
s_f0_fig = figure;

hold on
cmap = flipud(brewermap(length(mu_vec),'Spectral'));
colormap(cmap);

for m = 1:length(mu_vec)
    % draw plot indices
    f0_scatter_vec = 10.^s_f0_mu_titration(m).f0;
    s_scatter_vec = s_f0_mu_titration(m).sharpness/seq2;
    plot_options = find(s_scatter_vec>=0 & f0_scatter_vec >=1);
    plot_indices = randsample(plot_options,min([n_plot,length(plot_options)]),false);
    cw = mu_vec(m);
    cw = cw;
    % find the upper bound from the numerical scatters
    f0_axis = f0_vec(1:end-1) + diff(f0_vec); 
    s_num_vec = NaN(1,length(f0_vec)-1);    
    for f = 1:length(f0_vec)-1
        s_num_vec(f) = nanmax(s_scatter_vec(f0_scatter_vec < f0_vec(f+1)/alpha_factor & f0_scatter_vec >= f0_vec(f)/alpha_factor));
    end
%     s_bound_1 = (f0_vec./mu_vec(m)) ./ (1 + f0_vec./mu_vec(m)) .* seq2 .* (alpha_factor-1)./alpha_factor .* (alpha_factor + f0_vec) ./ (f0_vec- 1); 
    % NL: this bound is obtained in mathematica notebook entitled
    % "sharpness_vs_specificity_bound"
    s_bound = ((-1)+alpha_factor).^(-1).*(alpha_factor.^2+((-2)+alpha_factor).*f0_vec).*(cw+c.*f0_vec).^(-1).*seq2;    
%     s_bound_3 = ((-1)+alpha_factor).*alpha_factor.^(-1).*((-1)+f0_vec).^(-1).*f0_vec.*(alpha_factor+f0_vec).*(cw+c.*f0_vec).^(-1).* seq2;   
    
    scatter(f0_axis/alpha_factor, s_num_vec,...
          0.5*markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(m,:)); 
    
    plot(f0_vec/alpha_factor,s_bound/seq2,'--','Color',brighten(cmap(m,:),-0.5),'LineWidth',2)       
    
end

set(gca,'xscale','log')
h = colorbar;
h.Ticks = ([1 6 11]-0.5)/11;
h.TickLabels = {'10^{-2}','10^0','10^{2}'};
xlabel('fidelity (f^{neq}_0/f^{eq}_0)');
ylabel('effective sharpness (s^{neq}/s^{eq})')
% xlabel('activator fidelity (f_0/\alpha_factor)');
% ylabel('sharpness (dr/dc)')
ylabel(h, '\mu/\alpha_factor')


grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

s_f0_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

% saveas(s_f0_fig,[FigPath 's_vs_f0_vs_mu.png'])
% saveas(s_f0_fig,[FigPath 's_vs_f0_vs_mu.pdf'])

%% %% Track effective sharpness vs mu for f0 and s0-optimized networks %%%%

nBins = 51;
if ~exist([FigPath 's0_f0_motif_struct.mat'])    
    [s0_f0_motif_struct] = generate_s_vs_f0_vs_mu_scatters(nBins,sweep_options,alpha_factor,FigPath);
else
    load([FigPath 's0_f0_motif_struct.mat'],'s0_f0_motif_struct')
end

%%
seq = 1/4;
%%%%%%%%%%%%%%%%%
% generate predicted curves for the two motifs
mu_vec2 = s0_f0_motif_struct.mu_vec;
cw = 1./mu_vec2;

s_bound_fun2 = @(s0) alpha_factor.*cw.*s0.*(s0+((-1)+alpha_factor.*cw).*seq).*(s0+((-1)+2.*alpha_factor).*cw.* ...
                      s0+((-1)+cw.*(2+alpha_factor.*((-2)+alpha_factor.*cw))).*seq).^(-1);
                
s_bound_fun = @(f0) (f0./mu_vec2) ./ (1 + f0./mu_vec2) .* seq .* (alpha_factor-1) ./ alpha_factor .* (alpha_factor + f0) ./ (f0- 1); 
s_bound_s0 = s_bound_fun2(0.5);
s_bound_f0 = s_bound_fun(alpha_factor^2);

x_axis = (mu_vec2./alpha_factor) ;
%%%%%%%%%%%%%%%%%%%%%%%%
% make figure
close all
motif_fig = figure;
hold on
cmap = brewermap([],'Set2');


s1 = scatter(mu_vec2/alpha_factor, [s0_f0_motif_struct.s_max_vec],...
          0.5*markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor',cmap(9,:),'MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor','k'); 
% % sweep results for sharpness-optimized
% scatter(x_axis, [s0_f0_motif_struct.s_vec_spec],...
%           0.5*markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(4,:));
% prediction for sharpness-optimized
p1 = plot(x_axis,s_bound_f0,'--','Color',brighten(cmap(4,:),-.5),'LineWidth',3);

% % sweep for f0-optimized
% scatter(x_axis, [s0_f0_motif_struct.s_vec_sharp],...
%           0.5*markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha, 'MarkerFaceColor',cmap(5,:));        

% f0 prediction
p2 = plot(x_axis,s_bound_s0,'--','Color',brighten(cmap(5,:),-.5),'LineWidth',3);        

% optimally sharp network
% plot(mu_vec2/alpha_factor, [s_bound_s0(1:25) s_bound_f0(26:end)],'--','Color',[0 0 0 .2],'LineWidth',3);


set(gca,'xscale','log')
xlabel('\mu / \alpha_factor');
% set(gca, 'XDir','reverse')
ylabel('effective sharpness (s^{neq}/s^{eq})')
% ylabel(h, '\mu \alpha_factor')

legend([p2 p1 s1],'sharpness optimized','fidelity optimized','global optimum')

grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

motif_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(motif_fig,[FigPath 'motif_sharpness_plot.png'])
saveas(motif_fig,[FigPath 'motif_sharpness_plot.pdf'])



