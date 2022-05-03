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
p3Right_index = find(strcmp(metric_names,'p3Right'));
p3Wrong_index = find(strcmp(metric_names,'p3Wrong'));

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
tic
[sim_info_neq, sim_struct_neq] = param_sweep_multi_v2([sharpness_index cw_index],functionPath,sweep_options{:},...
                                          'half_max_flag',true,'equilibrium_flag',false,'specFactor',alpha_factor);
toc     

%%
tic
[sim_info_eq, sim_struct_eq] = param_sweep_multi_v2([sharpness_index cw_index],functionPath,sweep_options{:},...
                                          'half_max_flag',true,'equilibrium_flag',true,'specFactor',alpha_factor);
toc     

%%
metric_array = vertcat(sim_struct_neq.metric_array);
cwVec = metric_array(:,cw_index);
sharpnessVec = metric_array(:,sharpness_index);
pdRate = metric_array(:,rate_index);
p3RightProb = metric_array(:,p3Right_index);
p3WrongProb = metric_array(:,p3Wrong_index);
p4Prob = pdRate - p3RightProb - p3WrongProb;%metric_array(:,p3WrongIndex);

% find best networks for earch cw value and look at their state occupancy
% stats
cw_vec = linspace(nanmin(cwVec),nanmax(cwVec),51);
sharp_vec = NaN(1,length(cw_vec));
p3Right_vec = NaN(1,length(cw_vec));
p3Wrong_vec = NaN(1,length(cw_vec));
p4_vec = NaN(1,length(cw_vec));

for c = 2:length(cw_vec)
    % find index
    c_filter = cwVec<cw_vec(c)&cwVec>=cw_vec(c-1);
    [sharp_vec(c),si] = nanmax(sharpnessVec.*c_filter);
    % extract other values
    p3Right_vec(c) = p3RightProb(si);
    p3Wrong_vec(c) = p3WrongProb(si);
    p4_vec(c)= p4Prob(si);
end
    