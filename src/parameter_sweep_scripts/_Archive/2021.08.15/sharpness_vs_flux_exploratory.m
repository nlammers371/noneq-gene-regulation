% exploratory script to figure out weird disconnect between entropy rate,
% sharpness, and cycle time
clear 
close all
addpath(genpath('../utilities/'))


[~,metric_names] = calculateMetricsSym([]);
nStates = 4;
functionPath = ['../utilities/metricFunctions/n' num2str(nStates) '_OR/'];

% make sure we're linked to the appropriate function subfolder
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));

% define save path
OutPath = ['../../out/bivariate_parameter_sweeps_n' num2str(nStates) '_OR' filesep];
mkdir(OutPath);
                         

% get index of useful metrics
flux_index = find(strcmp(metric_names,'Flux'));
rate_index = find(strcmp(metric_names,'Production Rate'));
bin_index = find(strcmp(metric_names,'BinomialNoise'));
spec_alt_index = find(strcmp(metric_names,'specFactorAlt'));
sharp_right_index = find(strcmp(metric_names,'SharpnessRight'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
sharpness_norm_index = find(strcmp(metric_names,'SharpnessNormed'));
precision_index = find(strcmp(metric_names,'Precision'));
information_rate_index = find(strcmp(metric_names,'DecisionRateNorm'));
decision_time_index = find(strcmp(metric_names,'DecisionTimeNorm'));
phi_index = find(strcmp(metric_names,'Phi'));
affinity_index = find(strcmp(metric_names,'AffinityVec'));
dev_index = find(strcmp(metric_names,'deviationFactor'));
tau_index = find(strcmp(metric_names,'CycleTime'));


% set sim options
sweep_options = {'n_seeds',10,'n_iters_max',100,'nStates',nStates,'numCalcFlag',0};

% conduct sweep
tic;
[sim_info_neq, sim_struct_neq] = param_sweep_multi_v2([sharpness_index precision_index],functionPath,sweep_options{:},...
                                                            'half_max_flag',true,'equilibrium_flag',false,'saturationFlag',0);

[sim_info_eq, sim_struct_eq] = param_sweep_multi_v2([sharpness_index precision_index],functionPath,sweep_options{:},...
                                                            'half_max_flag',true,'equilibrium_flag',true,'saturationFlag',0);
% [sim_info_neq_slow, sim_struct_neq_slow] = param_sweep_multi_v2([sharpness_index precision_index],functionPath,sweep_options{:},...
%                                                             'half_max_flag',true,'equilibrium_flag',false,'TauCycleLimit',100);                                                          
toc;
%%

tau_norm = 5*60;

% make graph plot
metric_array_eq = vertcat(sim_struct_eq.metric_array);
rate_array_eq = vertcat(sim_struct_eq.rate_array);

metric_array = vertcat(sim_struct_neq.metric_array);
rate_array = vertcat(sim_struct_neq.rate_array);

phi_vec = metric_array(:,phi_index);
flux_vec = metric_array(:,flux_index);
tau_vec = metric_array(:,tau_index);
sharpness_vec = metric_array(:,sharpness_index);
sharpness_vec_eq = metric_array_eq(:,sharpness_index);
precision_vec = exp(metric_array(:,precision_index)).^2;

info_vec = metric_array(:,information_rate_index);
info_vec_eq = metric_array_eq(:,information_rate_index);

sharp_indices = find( sharpness_vec >=0.49);
precise_indices = find(precision_vec>=15.8 & sharpness_vec >=0);

inputParamsSharp = rate_array(sharp_indices(1),:);
inputParamsSharp(2:end) = inputParamsSharp(2:end) * tau_vec(sharp_indices(1))/tau_norm;
paramCellSharp = mat2cell(inputParamsSharp,size(inputParamsSharp,1),ones(1,size(inputParamsSharp,2)));
RNumSharp = RSymFun(paramCellSharp{:});

inputParamsPrecise = rate_array(precise_indices(1),:);
inputParamsPrecise(2:end) = inputParamsPrecise(2:end) * tau_vec(precise_indices(1))/tau_norm;
paramCellPrecise = mat2cell(inputParamsPrecise,size(inputParamsPrecise,1),ones(1,size(inputParamsPrecise,2)));
RNumPrecise = RSymFun(paramCellPrecise{:});


stateProbsSharp = steadyStateVecFunction(paramCellSharp{:}); 
% RTrueFluxNoisy = RNumSharp.*stateProbsNoisy;
% JNoisy = RTrueFluxNoisy(2,1) - RTrueFluxNoisy(1,2)

stateProbsPrecise = steadyStateVecFunction(paramCellPrecise{:}); 
% RTrueFluxPrecise = RNumPrecise.*stateProbsPrecise;
% JPrecise = RTrueFluxPrecise(2,1) - RTrueFluxPrecise(1,2)
%%
% sharp_rates = RNumSharp(RNumSharp>0);
% sharp_rates_norm = sharp_rates/sum(sharp_rates);
% precise_rates = RNumPrecise(RNumPrecise>0);
% precise_rates_norm = precise_rates / sum(precise_rates);

rate_array_norm = rate_array(:,2:end)./ sum(rate_array(:,2:end),2);
rate_array_norm_eq = rate_array_eq(:,2:end)./ sum(rate_array_eq(:,2:end),2);
paramCell = mat2cell(rate_array,size(rate_array,1),ones(1,size(rate_array,2)));
paramCellEq = mat2cell(rate_array_eq,size(rate_array_eq,1),ones(1,size(rate_array_eq,2)));
stateProbArray = steadyStateVecFunction(paramCell{:});
stateProbArrayEq = steadyStateVecFunction(paramCellEq{:});

% calculate rate and prob entropy measures
rate_entropy_vec = sum(rate_array_norm.*log(rate_array_norm),2);
rate_entropy_vec = rate_entropy_vec - min(rate_entropy_vec);
rate_entropy_vec = rate_entropy_vec/nanmax(rate_entropy_vec);

state_entropy_vec = sum(stateProbArray.*log(stateProbArray),2);
state_entropy_vec = state_entropy_vec - min(state_entropy_vec);
state_entropy_vec = state_entropy_vec/nanmax(state_entropy_vec);
% rate_entropy_vec_eq = -sum(rate_array_norm_eq.*log(rate_array_norm_eq),2);
% state_entropy_vec_eq = -sum(stateProbArrayEq.*log(stateProbArrayEq),2);
% sharp_state_entropy = stateProbsSharp*log(stateProbsSharp')
% precise_state_entropy = stateProbsPrecise*log(stateProbsPrecise')
% 
% sharp_rate_entropy = sharp_rates_norm'*log(sharp_rates_norm)
% precise_rate_entropy = precise_rates_norm'*log(precise_rates_norm)

% look for outliers
sharp_val = prctile(sharpness_vec,99.9);
precise_val = prctile(precision_vec,99.9);
info_val = prctile(info_vec,99.9);
sharp_filter = sharpness_vec >= sharp_val;
precise_filter = precision_vec >= precise_val;
info_filter = info_vec >= info_val;

figure;
hold on
scatter(state_entropy_vec,rate_entropy_vec)
% scatter(state_entropy_vec_eq,rate_entropy_vec_eq)
scatter(state_entropy_vec(sharp_filter),rate_entropy_vec(sharp_filter))
scatter(state_entropy_vec(precise_filter),rate_entropy_vec(precise_filter))
% scatter(state_entropy_vec(info_filter),rate_entropy_vec(info_filter))
%%
close all
gNoisy = digraph(RNumSharp>0);
gPrecise = digraph(RNumPrecise>0);

noisy_rates = RNumSharp(RNumSharp>0);
flux_rates = RNumSharp.*stateProbsNoisy;
flux_rates = flux_rates(flux_rates>0);

max_r = max(vertcat(noisy_rates,precise_rates));
min_r = min(vertcat(noisy_rates,precise_rates));
% noisy_rates = log10(noisy_rates) - log10(min_r) + 0.01;
% precise_rates = log10(precise_rates) - log10(min_r) + 0.01;

noisy_network = figure; 
hold on
% lw = log(RNumFlux(RNum>0)) - min(log(RNumFlux(RNum>0)))+1e-6;
p1 = plot(gNoisy,'LineWidth',noisy_rates);
p1.MarkerSize = 50*stateProbsNoisy;
p1.XData = [1 2 2 1];
p1.YData = [1 1 0 0];

flux = figure; 
hold on
% lw = log(RNumFlux(RNum>0)) - min(log(RNumFlux(RNum>0)))+1e-6;
p1 = plot(gNoisy,'LineWidth',flux_rates*50);
p1.MarkerSize = 50*stateProbsNoisy;
p1.XData = [1 2 2 1];
p1.YData = [1 1 0 0];


