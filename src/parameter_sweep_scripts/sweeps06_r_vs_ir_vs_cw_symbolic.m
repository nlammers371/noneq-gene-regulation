% script to conduct sweeps tracking IR vs pon for different CW values

% clear 
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

% Set Dropbox directory
DropboxFolder = 'S:\Nick\Dropbox\Nonequilibrium\Nick\SweepOutput';
writePath = [DropboxFolder filesep 'sweeps06_r_vs_info_vs_cw' filesep];
mkdir(writePath);


% get index of useful metrics
ir_index = find(strcmp(metric_names,'DecisionRateNorm'));
rate_index = find(strcmp(metric_names,'Production Rate'));
cw_index = find(strcmp(metric_names,'CW'));
rate_ent_index = find(strcmp(metric_names,'rateEntropy'));

% set sim options
sweep_options = {'n_seeds',5,'n_iters_max',50,'n_sim',5,'nStates',nStates};

% calculate sensitivity bound
alpha_factor = 100;

% specify plot options 
markerAlpha = 0.5; % marker transparency
markerSize = 75; % marker size
n_plot = 3e3;
bit_factor = log2(exp(1));

cw_vec = logspace(0,6,151);

%% Simulate perturbations to binding site of differing levels of severity
%%%%%%%%%%%%%%%%%%%%

try
    parpool(24)
catch
    % do nothing
end


parfor c = 1:length(cw_vec)
    % conduct sweeps    
    % global neq
    [sweep_info_neq, sweep_results_neq, sweep_results_neq_short] = param_sweep_multi_v3([rate_index ir_index],functionPath,sweep_options{:},...
                                                  'half_max_flag',false,'cw',cw_vec(c),...
                                                  'equilibrium_flag',false,'specFactor',alpha_factor,'paramBounds',paramBounds,'downsample_output',false);
                                            
    [sweep_info_eq, sweep_results_eq, sweep_results_eq_short] = param_sweep_multi_v3([rate_index ir_index],functionPath,sweep_options{:},...
                                                  'half_max_flag',false,'cw',cw_vec(c),...
                                                  'equilibrium_flag',true,'specFactor',alpha_factor,'paramBounds',paramBounds,'downsample_output',false);
                                                
%     [sweep_info_neq0, sweep_results_neq0] = param_sweep_multi_v3([rate_index rate_ent_index],functionPath,sweep_options{:},...
%                                                   'half_max_flag',false,'cw',cw_vec(c),...
%                                                   'equilibrium_flag',false,'specFactor',alpha_factor,'paramBounds',paramBounds);
                                                
    [sweep_info_eq0, sweep_results_eq0, sweep_results_eq0_short] = param_sweep_multi_v3([rate_index rate_ent_index],functionPath,sweep_options{:},...
                                                  'half_max_flag',false,'cw',cw_vec(c),...
                                                  'equilibrium_flag',true,'specFactor',alpha_factor,'paramBounds',paramBounds);
                                            
    % save
    suffix = ['cw' sprintf('%03d',c)];
    saveFun(suffix, writePath, sweep_info_neq, sweep_results_neq, sweep_info_eq, sweep_results_eq, sweep_info_eq0, sweep_results_eq0)    
end    

                                            