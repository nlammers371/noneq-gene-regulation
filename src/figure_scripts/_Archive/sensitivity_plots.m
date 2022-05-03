% Script to make plots for sensitivity figure

clear 
close all
addpath('../utilities/')

FigPath = '../../fig/sensitivity_plots/';
mkdir(FigPath);

% define physical constants
atp = 30.5 * 1000 / 6.022e23; % in joules
kT = 300*1.38e-23;

green = [191 213 151]/256;
% (1) plot sensitivity vs. gamma at equilibrium 
gamma_vec = linspace(0,1*atp);
sensitivity_vec = 0.25 * (1-exp(-gamma_vec./kT)) ./ (1+exp(-gamma_vec./kT));

gamma_fig = figure;
hold on
area(gamma_vec/atp, sensitivity_vec,'FaceColor',green,'EdgeColor',green,'LineWidth',1);
plot(gamma_vec/atp, sensitivity_vec,'-','Color','k','LineWidth',1.5);
p = plot(0,0);
StandardFigurePBoC(p,gca)

xlabel('\gamma (ATP units)')
ylabel('sensitivity (dr/dc)')
set(gca,'FontSize',14)

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
ylim([0 0.3])
gamma_fig.InvertHardcopy = 'off';
set(gca,'Layer','top')

saveas(gamma_fig,[FigPath 'sensitivity_vs_gamma.png'])
saveas(gamma_fig,[FigPath 'sensitivity_vs_gamma.pdf'])

% calculate implied occupancy of states 1 and 3
occ = exp(-0.5*atp/kT) / (2*exp(-0.5*atp/kT) + 2)