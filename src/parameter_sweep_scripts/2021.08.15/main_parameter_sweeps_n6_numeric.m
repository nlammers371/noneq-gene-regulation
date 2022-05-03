% script to call core parameter sweep function to examine tradeoffs between
% different network behaviors

clear 
% close all
addpath(genpath('../utilities/'))


[~,metric_names] = calculateMetricsNumeric([]);
nStates = 6;
functionPath = ['../utilities/metricFunctions/n' num2str(nStates) '_OR_NUM/'];

% make sure we're linked to the appropriate function subfolder
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));

% define save path
OutPath = ['../../out/bivariate_parameter_sweeps_n' num2str(nStates) '_numeric' filesep];
mkdir(OutPath);
                         

% get index of useful metrics
flux_index = find(strcmp(metric_names,'Flux'));
rate_index = find(strcmp(metric_names,'Production Rate'));
spec_index = find(strcmp(metric_names,'Specificity'));
spec_alt_index = find(strcmp(metric_names,'specFactorAlt'));
sharp_right_index = find(strcmp(metric_names,'SharpnessRight'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
sharpness_norm_index = find(strcmp(metric_names,'SharpnessNormed'));
precision_index = find(strcmp(metric_names,'Precision'));
decision_rate_index = find(strcmp(metric_names,'DecisionRateNorm'));
decision_time_index = find(strcmp(metric_names,'DecisionTimeNorm'));
phi_index = find(strcmp(metric_names,'Phi'));
affinity_index = find(strcmp(metric_names,'AffinityVec'));
dev_index = find(strcmp(metric_names,'deviationFactor'));
cw_index = find(strcmp(metric_names,'CW'));


% set sim options
sweep_options = {'n_sim',1,'n_seeds',5,'n_iters_max',50,'nStates',nStates,'numCalcFlag',1};

%% %%%%%%%%%%%%%%%% info rate vs energy flux per cycle %%%%%%%%%%%%%%%%%%%%
tic
% [simInfoEq, sim_struct_eq] = param_sweep_multi_v2([spec_index decision_rate_index],functionPath, sweep_options{:},...
%                                           'half_max_flag',true,'equilibrium_flag',true,'cw',1);
                                        
[simInfoNeq, sim_struct_neq] = param_sweep_multi_v2([cw_index sharpness_norm_index],functionPath, sweep_options{:},...
                                          'half_max_flag',true,'equilibrium_flag',false,'cw',1);                                        
toc     

%%
close all
[~,mi] = nanmax(sim_struct_neq(1).metric_array(:,decision_rate_index));

inputParams = sim_struct_neq(1).rate_array(mi,:);
paramCell = mat2cell(inputParams,size(inputParams,1),ones(1,size(inputParams,2)));
RNum = RSymFun(paramCell{:});

numerical_precision = 10;
state_probs = calculate_ss_num(RNum,numerical_precision); 
RNumFlux = RNum.*state_probs';
RNumFlux = RNumFlux / max(RNumFlux(:)) * 10;

RNumGraph = RNum;
% RNumGraph(RNumFlux<=1e-4) = 0;

g1 = digraph(RNumGraph>0);

RNumGraphSmall = RNum; 
RNumGraphSmall (RNumFlux>1e-4) = 0;

g2 = digraph(RNumGraphSmall>0);

figure; 

hold on
% lw = log(RNumFlux(RNum>0)) - min(log(RNumFlux(RNum>0)))+1e-6;
p1 = plot(g1,'LineWidth',RNumFlux(RNum>0));
p1.MarkerSize = 50*state_probs;
p1.XData = [1 2 2 1 0 0];
p1.YData = [1 1 0 0 0 1 ];


cw_vec = sim_struct_neq(1).metric_array(:,cw_index);
[~,mi2] = nanmax(sim_struct_neq(1).metric_array(cw_vec>=3.85,decision_rate_index));

inputParams2 = sim_struct_neq(1).rate_array(mi2,:);
paramCell2 = mat2cell(inputParams2,size(inputParams2,1),ones(1,size(inputParams2,2)));
RNum2 = RSymFun(paramCell2{:});

state_probs2 = calculate_ss_num(RNum2,numerical_precision); 
RNumFlux2 = RNum2;%.*state_probs2';
RNumFlux2 = RNumFlux2 / 1e4 * 25;

g1 = digraph(RNum>0);

figure; 

hold on
% lw = log(RNumFlux(RNum>0)) - min(log(RNumFlux(RNum>0)))+1e-6;
p1 = plot(g1,'LineWidth',RNumFlux2(RNum2>0));
p1.MarkerSize = 50*state_probs2;
p1.XData = [1 2 2 1 0 0];
p1.YData = [1 1 0 0 0 1 ];

% p2 = plot(g2,'LineWidth',RNumFlux(RNumGraphSmall>0));
% p2.LineStyle = '--';
% p2.MarkerSize = 50*state_probs;
% p2.XData = [1 2 2 1 0 0];
% p2.YData = [1 1 0 0 0 1 ];
