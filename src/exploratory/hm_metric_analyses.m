% script to track observed sharpness as a function of CW
clear 
close all
addpath(genpath('../utilities/'))

% set basic parameters
nStates = 6;
paramBounds = repmat([-10 ; 6],1,11); % constrain transition rate magnitude
[~,~,metric_names] = calculateMetricsSym_v2([]);

% specify function path
functionPath = ['../utilities/metricFunctions/n' num2str(nStates) '_OR/'];

% make sure we're linked to the appropriate function subfolder% make sure we're linked to the appropriate function subfolder
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));

% define save path
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\manuscript\';
% DropboxFolder = 'S:\Nick\Dropbox\Nonequilibrium\Nick\manuscript\';

FigPath = [DropboxFolder 'exploratory_analyese' filesep];
mkdir(FigPath);         

% get index of useful metrics
spec_index = find(strcmp(metric_names,'Specificity'));
sharp_right_norm_index = find(strcmp(metric_names,'SharpnessRightNorm'));
ir_index = find(strcmp(metric_names,'DecisionRateNorm'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
precision_index = find(strcmp(metric_names,'Precision'));
rate_index = find(strcmp(metric_names,'Production Rate'));
precision_right_index = find(strcmp(metric_names,'PrecisionRight'));
cw_index = find(strcmp(metric_names,'CW'));
tau_index = find(strcmp(metric_names,'CycleTime'));

% set sim options
sweep_options = {'n_seeds',5,'n_iters_max',50,'n_sim',1,'nStates',nStates};

% calculate sensitivity bound
alpha_factor = 100;
% seq = 1/4;

% specify plot options 
markerAlpha = 0.5; % marker transparency
markerSize = 75; % marker size
% n_plot = 3e3;
% bit_factor = log2(exp(1));
%% 
% observed sharpness vs cw
tic
[sim_info_neq, sim_struct_neq] = param_sweep_multi_v3([rate_index sharpness_index],functionPath,sweep_options{:},...
                                          'half_max_flag',false,...
                                          'equilibrium_flag',false,'cw',1e4,'specFactor',alpha_factor,'paramBounds',paramBounds);
                                       
toc     

%% Let's pull out a random set of rates and see what happens when we fiddle with stuff

s_filter = sim_struct_neq(1).metric_array(:,sharpness_index)>=.05 & sim_struct_neq(1).metric_array(:,rate_index)<=.25;
sharp_vec_raw = sim_struct_neq(1).metric_array(s_filter,sharpness_index);
prec_vec_raw = exp(sim_struct_neq(1).metric_array(s_filter,precision_index));
rate_vec_raw = sim_struct_neq(1).metric_array(s_filter,rate_index);
tau_vec_raw = sim_struct_neq(1).metric_array(s_filter,tau_index);

sharp_hm_pd = sharp_vec_raw ./(rate_vec_raw.*(1-rate_vec_raw))/4;
prec_hm_pd = prec_vec_raw .*(rate_vec_raw.*(1-rate_vec_raw)).^2 *16;

sharp_hm_true = NaN(size(sharp_hm_pd));
sharp_max_true = NaN(size(sharp_hm_pd));
precision_hm_true = NaN(size(sharp_hm_pd));

rate_array = sim_struct_neq(1).rate_array(s_filter,:);

% now calculate induction curves
c_vec = logspace(-3,3,1e3);
induction_curve_array = NaN(length(c_vec),sum(s_filter));

for i = 1:length(sharp_hm_true)    

    neq_rates = repmat(rate_array(i,:),length(c_vec),1);   
    neq_rates(:,1) = c_vec';

    % convert to cell array and run calculations
    paramCellNeq = mat2cell(neq_rates,size(neq_rates,1),ones(1,size(neq_rates,2)));

    % production rate 
    r_curve = productionRateFunction(paramCellNeq{:});
    induction_curve_array(:,i) = r_curve;
    p_curve = intrinsicVarianceFunction(paramCellNeq{:});

    TauOn = TauONFunction(paramCellNeq{:});
    TauOff = TauOFFFunction(paramCellNeq{:});
    TauCycle = TauOff+TauOn;   
    % find HM point and calculate sharpness
    [~,hm_i] = min(abs(r_curve-0.5));
    if ~isempty(hm_i) && nanmax(r_curve)>0.98 && nanmin(r_curve)<0.02
        sharp_vec = diff(r_curve')./diff(c_vec) .* (c_vec(1:end-1)+diff(c_vec)/2);
        sharp_hm_true(i) = sharp_vec(hm_i-1);
        precision_hm_true(i) = TauCycle(hm_i)./p_curve(hm_i);
        sharp_max_true(i) = nanmax(sharp_vec);
    end
end

