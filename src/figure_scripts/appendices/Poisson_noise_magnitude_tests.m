% This script performs simple calculations to examine the expected
% contribution from Poisson nosie in mRNA production to the overall amount
% of observed transcriptional noise. 
clear
close all
addpath(genpath('../utilities'))
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\';
FigPath = [DropboxFolder '\manuscript\appendices' filesep];
mkdir(FigPath)
%%
close all
PolII_initiation_rate = 20; % Pol II per minute (based off of estimates for eve stripe 2 from Lammers 2020)
% leverage fact that variance per burst cycle is equal to some constant times b^2, where
% b=a(a-1). For the case of the 4 state system away from equalibrium, this
% constant is just 1, such that:
%      \sigma <= r^2*\tau_b x b^2 x 1
% The key is to realize that the relative contribution from poisson remains
% constant with a, while bursting component increases linearly with tau_b

noise_prefactor = 1; % 4 state system at equilibrium
burst_time_vec = 1:1:(6*60); % burst cycle times (in minutes)
a_vec = linspace(0.02,0.98)';

poisson_noise_array = repmat(a_vec.*PolII_initiation_rate,1,length(burst_time_vec));
burst_noise_array = PolII_initiation_rate^2*(a_vec.*(1-a_vec)).^2 .* burst_time_vec .* noise_prefactor;
rel_array = poisson_noise_array./(poisson_noise_array+burst_noise_array);

poisson_noise_fig = figure;
cmap = brewermap([],'YlOrRd');
colormap(cmap)
p = pcolor(poisson_noise_array);
p.EdgeAlpha = 0.05;

set(gca,'Fontsize',14);
xlabel('burst cycle time (minutes)')
ylabel('fraction of time in active state (a)')

h = colorbar;
ylabel(h,'predicted Poisson noise component (\sigma_{mRNA}^2)')
yt = yticks;
set(gca,'yticklabels',round(a_vec(yt),2)*100)
xt = xticks;
set(gca,'xticklabels',round(burst_time_vec(xt)/5)*5)
poisson_noise_fig.Renderer='Painters';
saveas(poisson_noise_fig,[FigPath 'poisson_heatmap.png']);
saveas(poisson_noise_fig,[FigPath 'poisson_heatmap.pdf']);

% burst noise
burst_noise_fig = figure;
colormap(cmap)
p = pcolor(burst_noise_array);
p.EdgeAlpha = 0.05;

set(gca,'Fontsize',14);
xlabel('burst cycle time (minutes)')
ylabel('fraction of time in active state (a)')

h = colorbar;
ylabel(h,'minimum bursting noise component (\sigma_{burst}^2)')
yt = yticks;
xt = xticks;
set(gca,'xticklabels',round(burst_time_vec(xt)/5)*5)
burst_noise_fig.Renderer='Painters';
set(gca,'yticklabels',round(a_vec(yt),2)*100)
saveas(burst_noise_fig,[FigPath 'burst_noise_heatmap.png']);
saveas(burst_noise_fig,[FigPath 'burst_noise_heatmap.pdf']);

% relative poisson contribution
rel_noise_fig = figure;
cmap = flipud(brewermap([],'Spectral'));
colormap(cmap)
p = pcolor(rel_array);
p.EdgeAlpha = 0.05;

set(gca,'Fontsize',14);
xlabel('burst cycle time (minutes)')
ylabel('fraction of time in active state (a)')
caxis([0 1])
h = colorbar;
rel_noise_fig.Renderer='Painters';
ylabel(h,'\sigma_{burst}^2/(\sigma_{mRNA}+\sigma_{burst})')
yt = yticks;
xt = xticks;
set(gca,'xticklabels',round(burst_time_vec(xt)/5)*5)
set(gca,'yticklabels',round(a_vec(yt),2)*100)
saveas(rel_noise_fig,[FigPath 'rel_noise_heatmap.png']);
saveas(rel_noise_fig,[FigPath 'rel_noise_heatmap.pdf']);

