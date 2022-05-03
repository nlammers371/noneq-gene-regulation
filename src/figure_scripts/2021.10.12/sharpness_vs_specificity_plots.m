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

FigPath = [DropboxFolder 'sharpness_vs_specificity' filesep];
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
cw = 1e-6; % precise value unimportant...just something small enough to be negligible
tic
[sim_info_neq, sim_struct_neq] = param_sweep_multi_v2([sharpness_index spec_index],functionPath,sweep_options{:},...
                                          'half_max_flag',true,'cw',cw,...
                                          'equilibrium_flag',false,'specFactor',alpha_factor);

[sim_info_eq, sim_struct_eq] = param_sweep_multi_v2([sharpness_index spec_index],functionPath,sweep_options{:},...
                                          'half_max_flag',true,'cw',cw,...
                                          'equilibrium_flag',true,'specFactor',alpha_factor);                                        
toc     

close all

%% cw = 1;
% get predicted tradeoff bound
% s0_bound = seq+(alpha_factor.^2+(-1).*f0_vec).*cw.*((-1).*alpha_factor+f0_vec+((-1)+alpha_factor).*f0_vec.*cw).^( ...
%   -1).*seq;
seq = 1/4;
cr = 1;
n_plot1 = 1e4;
s0_bound = (cr.^(-1)+((-1)+alpha_factor).^(-1).*cr.^(-1).*(alpha_factor.^2+(-1).*f0_vec).*f0_vec.^(-1));             
s0_bound2 = (cr.^(-1)+(alpha_factor.^2+(-1).*f0_vec).*((-1).*alpha_factor.*cw+(((-1)+alpha_factor).*cr+cw).*f0_vec).^(-1));

% generate vectors to plot (neq)
results_array_neq = vertcat(sim_struct_neq.metric_array);
f0_scatter_vec_neq = 10.^results_array_neq(:,spec_index);
s0_scatter_vec_neq = results_array_neq(:,sharpness_index)/seq;
plot_options_neq = find(s0_scatter_vec_neq>=0 & f0_scatter_vec_neq >= 1);
plot_indices_neq = randsample(plot_options_neq,min([n_plot1,length(plot_options_neq)]),false);
% generate vectors to plot (eq)
results_array_eq = vertcat(sim_struct_eq.metric_array);
f0_scatter_vec_eq = 10.^results_array_eq(:,spec_index);
s0_scatter_vec_eq = results_array_eq(:,sharpness_index)/seq;
plot_options_eq = find(s0_scatter_vec_eq>=0 & f0_scatter_vec_eq >= 1);
plot_indices_eq = randsample(plot_options_eq,min([n_plot1,length(plot_options_eq)]),false);

% make figure
s0_f0_fig = figure;
hold on
cmap = brewermap([],'Set2');

% sneq = scatter(f0_scatter_vec_neq(plot_indices_neq)*alpha_factor,s0_scatter_vec_neq(plot_indices_neq),...
%       markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha*0.5, 'MarkerFaceColor',cmap(2,:));
    
seq = scatter(f0_scatter_vec_eq(plot_indices_eq)*alpha_factor,s0_scatter_vec_eq(plot_indices_eq),...
      markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha*0.5, 'MarkerFaceColor',cmap(3,:));    
    
p = plot(f0_vec,s0_bound2,'-.','Color','k','LineWidth',3);
ylim([0 2])
xlim([alpha_factor alpha_factor^2])
set(gca,'xscale','log')    
xlabel('intrinsic specificity (f_0)');
ylabel('intrinsic sharpness (s_0/s_0^{eq})')
% grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
set(gca,'xtick',[alpha_factor alpha_factor^1.5 alpha_factor^2],'xticklabels',{'\alpha','\alpha^{1.5}','\alpha^2'})
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

s0_f0_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
set(p,'Visible','off')

saveas(s0_f0_fig,[FigPath 's0_vs_f0.png'])

% now add non-equilibrium results
sneq = scatter(f0_scatter_vec_neq(plot_indices_neq)*alpha_factor,s0_scatter_vec_neq(plot_indices_neq),...
      markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha*0.5, 'MarkerFaceColor',cmap(2,:));
uistack(seq,'top')
saveas(s0_f0_fig,[FigPath 's0_vs_f0_neq.png'])

set(p,'Visible','on')
uistack(p,'top')
saveas(s0_f0_fig,[FigPath 's0_vs_f0_neq_bound.png'])

delete(seq)
delete(sneq)
set(p,'Visible','on')
grid on
saveas(s0_f0_fig,[FigPath 's0_vs_f0_neq.pdf'])
saveas(s0_f0_fig,[FigPath 's0_vs_f0_neq.pdf'])



%% %%%%%%%%%%%%%%%% sharpness vs activator fidelity %%%%%%%%%%%%%%%%%%%%%
close all

% [sim_info_neq, sim_struct_neq] = param_sweep_multi([sharpness_index cw_index],sweep_options{:},...
%                                                   'half_max_flag',true,'wrongFactorConcentration',mu_vec(m),'equilibrium_flag',false);

% calculate sensitivity bound
alpha_factor = sim_info_neq.specFactor;
f0_vec = logspace(log10(alpha_factor),log10(alpha_factor^2),101);
mu_vec = [1 10 100 1e3 1e4];%logspace(log10(1),log10(alpha_factor^2),11);


sharpness_cell = cell(1,length(mu_vec));
s0_cell = cell(1,length(mu_vec));
spec_vec = cell(1,length(mu_vec));

s_f0_mu_titration = struct;

if ~exist([FigPath 's_f0_mu_titration.mat'])
    for  m = 1:length(mu_vec)
        tic
        [sim_info_neq, sim_struct_neq] = param_sweep_multi_v2([sharpness_index spec_index],functionPath, sweep_options{:},...
                                                  'half_max_flag',true,'cw',mu_vec(m),'equilibrium_flag',false);

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
cr = 1;
n_plot = 1e5;
% make figure
s_f0_fig = figure;

hold on
cmap = flipud(brewermap(length(mu_vec),'Spectral'));
colormap(cmap);

set(gca,'xscale','log','xtick',[1 10 100],'xticklabels',{'\alpha','\alpha^{1.5}','\alpha^2'})
h = colorbar;
h.Ticks = 0.1:.2:1;
% h.TickLabels = {'10^{-2}','10^{-1}','10^0','10^{1}','10^{2}'};
h.TickLabels = {'\alpha^0','\alpha^{0.5}','\alpha','\alpha^{1.5}','\alpha^2'};
xlabel('intrinsic specificity (f_0)');
ylabel('effective sharpness (s_*)')
% xlabel('activator fidelity (f_0/\alpha_factor)');
% ylabel('sharpness (dr/dc)')
ylabel(h, 'c_w/c_r')

% grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

s_f0_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

ylim([0 0.5])
xlim([1e0 1e2])

s_extra = {};
s_keep = {};
p_keep = {};
for m = 1:length(mu_vec)
    % draw plot indices
    f0_scatter_vec_neq = 10.^s_f0_mu_titration(m).f0;
    s_scatter_vec = s_f0_mu_titration(m).sharpness;
    plot_options_neq = find(s_scatter_vec>=0 & f0_scatter_vec_neq >=1);
    plot_indices_neq = randsample(plot_options_neq,min([n_plot,length(plot_options_neq)]),false);
    cw = mu_vec(m);
%     cw = cw;
    % find the upper bound from the numerical scatters
    f0_axis = f0_vec(1:end-1) + diff(f0_vec); 
    s_num_vec = NaN(1,length(f0_vec)-1);    
    for f = 1:length(f0_vec)-1
        s_num_vec(f) = nanmax(s_scatter_vec(f0_scatter_vec_neq < f0_vec(f+1)/alpha_factor & f0_scatter_vec_neq >= f0_vec(f)/alpha_factor));
    end
%     s_bound_1 = (f0_vec./mu_vec(m)) ./ (1 + f0_vec./mu_vec(m)) .* seq2 .* (alpha_factor-1)./alpha_factor .* (alpha_factor + f0_vec) ./ (f0_vec- 1); 
    % NL: this bound is obtained in mathematica notebook entitled
    % "sharpness_vs_specificity_bound"
    s_bound = ((-1)+alpha_factor).^(-1).*(alpha_factor.^2+((-2)+alpha_factor).*f0_vec).*(cw+cr.*f0_vec).^(-1).*seq2;    
    
    s_extra{end+1} = scatter(f0_scatter_vec_neq(plot_indices_neq), s_scatter_vec(plot_indices_neq),...
                    0.75*markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',.75, 'MarkerFaceColor',brighten(cmap(m,:),0.25));        
        
    s_keep{end+1} = scatter(f0_axis(1:2:length(f0_axis))/alpha_factor, s_num_vec(1:2:length(f0_axis)),...
                    0.75*markerSize,'MarkerEdgeAlpha',.5,'MarkerEdgeColor','k','MarkerFaceAlpha',0.75, 'MarkerFaceColor',brighten(cmap(m,:),0.25));
                
    if m ~= 3
        p_keep{end+1} = plot(f0_vec/alpha_factor,s_bound,'-.','Color',brighten(cmap(m,:),-0.5),'LineWidth',3);       
    else
        p_keep{end+1} = plot(f0_vec/alpha_factor,s_bound,'-.','Color',brighten(cmap(m,:),-0.8),'LineWidth',3);
    end
    
    saveas(s_f0_fig,[FigPath 's_vs_f0_vs_cw_scatter_' sprintf('%02d',m) '.png'])
end

for m = 1:length(mu_vec)
    set(s_keep{m},'Visible','off')
    set(p_keep{m},'Visible','off')
end

saveas(s_f0_fig,[FigPath 's_vs_f0_vs_cw_scatter.png'])

% Make second version without background scatters
for m = 1:length(mu_vec)
%     set(s_keep{m},'Visible','on')
    set(p_keep{m},'Visible','on')
    delete(s_extra{m})
end
grid on
saveas(s_f0_fig,[FigPath 's_vs_f0_vs_cw.png'])
saveas(s_f0_fig,[FigPath 's_vs_f0_vs_cw.pdf'])

%% %% Track effective sharpness vs mu for f0 and s0-optimized networks %%%%

tic
[sim_info_cw_neq, sim_struct_cw_neq] = param_sweep_multi_v2([sharpness_index cw_index],functionPath,sweep_options{:},...
                                          'half_max_flag',true,'equilibrium_flag',false,'specFactor',alpha_factor);

[sim_info_cw_eq, sim_struct_cw_eq] = param_sweep_multi_v2([sharpness_index cw_index],functionPath,sweep_options{:},...
                                          'half_max_flag',true,'equilibrium_flag',true,'specFactor',alpha_factor);                                        
toc     

%%
seq = 1/4;
cr = 1;
n_plot3 = 1e3;
% extract vectors to plot (neq)
results_array_neq = vertcat(sim_struct_cw_neq.metric_array);
cw_scatter_vec_neq = results_array_neq(:,cw_index);
s_scatter_vec_neq = results_array_neq(:,sharpness_index)/seq;
plot_options_cw_neq = find(s_scatter_vec_neq>=0);
plot_indices_cw_neq = randsample(plot_options_cw_neq,min([n_plot,length(plot_options_cw_neq)]),false);

% extract vectors to plot (eq)
results_array_eq = vertcat(sim_struct_cw_eq.metric_array);
cw_scatter_vec_eq = results_array_eq(:,cw_index);
s_scatter_vec_eq = results_array_eq(:,sharpness_index)/seq;
plot_options_cw_eq = find(s_scatter_vec_eq>=0);
plot_indices_cw_eq = randsample(plot_options_cw_eq,min([n_plot,length(plot_options_cw_eq)]),false);


%%%%%%%%%%%%%%%%%
% generate predicted curves for the two motifs
cw_vec = logspace(nanmin(cw_scatter_vec_neq),nanmax(cw_scatter_vec_neq),101);
% cw = 1./cw_vec;
s_bound_fun = @(f0) ((-1)+alpha_factor).^(-1).*(alpha_factor.^2+((-2)+alpha_factor).*f0).*(cw_vec+cr.*f0).^(-1);
                
% s_bound_fun = @(f0) (f0./cw_vec) ./ (1 + f0./cw_vec) .* seq .* (alpha_factor-1) ./ alpha_factor .* (alpha_factor + f0) ./ (f0- 1); 
s_bound_s0 = s_bound_fun(alpha_factor);
s_bound_f0 = s_bound_fun(alpha_factor^2);
s_bound_eq = alpha_factor ./ (alpha_factor + cw_vec);
cw_axis = (cw_vec) ;

%%%%%%%%%%%%%%%%%%%%%%%%
% make figure
close all
motif_fig = figure;
hold on
cmap = brewermap([],'Set2');


% plot numerical results
sneq = scatter(10.^cw_scatter_vec_neq(plot_indices_cw_neq),s_scatter_vec_neq(plot_indices_cw_neq),...
      markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha*0.5, 'MarkerFaceColor',cmap(2,:));
    
seq = scatter(10.^cw_scatter_vec_eq(plot_indices_cw_eq),s_scatter_vec_eq(plot_indices_cw_eq),...
      markerSize,'MarkerEdgeAlpha',.25,'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha*0.5, 'MarkerFaceColor',cmap(3,:));    


% s0-optimized prediction
p1 = plot(cw_axis,s_bound_f0,'-.','Color',brighten(cmap(4,:),-.75),'LineWidth',3);
% f0-optimized prediction
p2 = plot(cw_axis,s_bound_s0,'-.','Color',brighten(cmap(5,:),-.75),'LineWidth',3);   
% equilibrium bound
p3 = plot(cw_axis,s_bound_eq,'-.','Color',brighten(cmap(3,:),-.75),'LineWidth',3);   


set(gca,'xscale','log')
xlabel('c_w / c_r');
ylabel('effective sharpness (s_*/s_0^{eq})')

set(gca,'xtick',[1 1e2 1e4 1e6],'xticklabels',{'\alpha^0','\alpha^1','\alpha^2','\alpha^3'})
% legend([p2 p1 s1],'sharpness optimized','fidelity optimized','global optimum')

xlim([1e0 1e4])
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

motif_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

set(p3,'Visible','off')
saveas(motif_fig,[FigPath 'effective_sharpness_vs_cw_scatter_s0_f0.png'])
set(p2,'Visible','off')
saveas(motif_fig,[FigPath 'effective_sharpness_vs_cw_scatter_s0.png'])
set(p1,'Visible','off')


saveas(motif_fig,[FigPath 'effective_sharpness_vs_cw_scatter.png'])

delete(seq)
delete(sneq)

set(p1,'Visible','on')
set(p2,'Visible','on')
set(p3,'Visible','on')

grid on
saveas(motif_fig,[FigPath 'effective_sharpness_vs_cw.png'])
saveas(motif_fig,[FigPath 'effective_sharpness_vs_cw.pdf'])
% saveas(motif_fig,[FigPath 'motif_sharpness_plot.png'])
% saveas(motif_fig,[FigPath 'motif_sharpness_plot.pdf'])



