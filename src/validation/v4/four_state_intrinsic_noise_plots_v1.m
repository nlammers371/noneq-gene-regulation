% This script plots results of intrinsic noise simulations agains
% theoretical expectations

clear
close all
addpath('../../utilities')

% define paths
ReadPath = '../../../out/validation_v4/';
mkdir(ReadPath);
FigPath = '../../../fig/validation_v4/';
mkdir(FigPath);

% specify project to read in
sim_label = 'intrinsic_sim_v1';

% load results
load([ReadPath sim_label '_simulation_results.mat'])

% iterate through rate array and calculate theoretical expectation in each
% case

% sim info
rate_array = sim_struct.rate_array;
c_vec = sim_struct.c_vec;

% sim results
var_vec_sim = sim_struct.var_r_vec;
pd_rate_vec_sim = sim_struct.mean_r_vec;

% initialize vectors
var_vec_theory = NaN(size(c_vec));
pd_rate_vec_theory = NaN(size(c_vec));
pd_rate_vec_theory_alt = NaN(size(c_vec));
           
for i = 1:length(c_vec)
  % calculate production rate and sharpness     
  pd_rate_vec_theory(i) = fourStateProductionGeneral(rate_array(i,:),c_vec(i));      
  var_vec_theory(i) = fourStateVarianceGeneral(rate_array(i,:),c_vec(i));
end
  
% make figures
mSize = 50;
mean_fig = figure;
cmap = brewermap(9,'Set2');
scatter(pd_rate_vec_theory,pd_rate_vec_sim,mSize,'MarkerEdgeColor','k','MarkerFaceAlpha',.5, 'MarkerFaceColor',cmap(2,:))
grid on
box on    
xlabel('theoretical prediction (r)')
ylabel('simulation result (r)')
set(gca,'FontSize',14)
% set(gca,'yScale','log')
% set(gca,'xScale','log')
% legend('non-equilibrium','equilibrium', 'Location','northeast') 
saveas(mean_fig,[FigPath 'mean_rate_validation.png'])
saveas(mean_fig,[FigPath 'mean_rate_validation.pdf'])

var_fig = figure;
scatter(var_vec_theory,var_vec_sim,mSize,'MarkerEdgeColor','k','MarkerFaceAlpha',.5, 'MarkerFaceColor',cmap(3,:))
grid on
box on    
xlabel('theoretical prediction (\sigma_{int})')
ylabel('simulation result (\sigma_{int})')
set(gca,'FontSize',14)
set(gca,'yScale','log')
set(gca,'xScale','log')
% legend('non-equilibrium','equilibrium', 'Location','northeast') 
saveas(var_fig,[FigPath 'intrinsic_noise_validation.png'])
saveas(var_fig,[FigPath 'intrinsic_noise_validation.pdf'])