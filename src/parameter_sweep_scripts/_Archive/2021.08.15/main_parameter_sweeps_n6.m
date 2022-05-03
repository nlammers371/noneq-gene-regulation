% script to call core parameter sweep function to examine tradeoffs between
% different network behaviors

clear 
close all
addpath(genpath('../utilities/'))

% set basic parameters
nStates = 6;
rate_bounds = repmat([-6 ; 6],1,3*nStates-4); % constrain transition rate magnitude
[~,metric_names] = calculateMetricsSym([]);

% make sure we're linked to the appropriate function subfolder
functionPath = ['../utilities/metricFunctions/n' num2str(nStates) '_OR/'];
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));

% define save path
OutPath = ['../../out/bivariate_parameter_sweeps_n' num2str(nStates) filesep];
mkdir(OutPath);
                         

% get index of useful metrics
flux_index = find(strcmp(metric_names,'Flux'));
rate_index = find(strcmp(metric_names,'Production Rate'));
spec_index = find(strcmp(metric_names,'Specificity'));
spec_alt_index = find(strcmp(metric_names,'specFactorAlt'));
sharp_right_index = find(strcmp(metric_names,'SharpnessRight'));
sharp_right_norm_index = find(strcmp(metric_names,'SharpnessRightNorm'));
sharp_spec_index = find(strcmp(metric_names,'SharpnessFactor'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
sharpness_norm_index = find(strcmp(metric_names,'SharpnessNormed'));


precision_index = find(strcmp(metric_names,'Precision'));
decision_rate_index = find(strcmp(metric_names,'DecisionRateNorm'));
phi_index = find(strcmp(metric_names,'Phi'));
decision_time_index = find(strcmp(metric_names,'DecisionTimeNorm'));
affinity_index = find(strcmp(metric_names,'AffinityVec'));
dev_index = find(strcmp(metric_names,'deviationFactor'));
cw_index = find(strcmp(metric_names,'CW'));


% set sim options
sweep_options = {'n_seeds',10,'n_iters_max',50,'nStates',nStates,'numCalcFlag',0};

%% %%%%%%%%%%%%%%%% info rate vs energy flux per cycle %%%%%%%%%%%%%%%%%%%%
tic
% [sim_info_eq, sim_struct_eq] = param_sweep_multi_v2([sharpness_index dev_index],functionPath,sweep_options{:},...
%                                           'half_max_flag',false,'cw',1,'equilibrium_flag',true);
cw = 1e4;                                        
[simInfoNeq, sim_struct_neq] = param_sweep_multi_v2([spec_index decision_rate_index],functionPath,sweep_options{:},...
                                          'half_max_flag',true,'cw',cw,'equilibrium_flag',false);                                        
toc     
%%
seq = 0.25;
mo = cw;
beta = 100;
f0 = logspace(log10(beta),log10(beta^2));
sneq = (seq + ((beta^2 - f0) *seq)./(-beta * mo + f0 *(-1 + beta + mo)))/seq;


%%
cw_vec = sim_struct_neq(1).metric_array(:,cw_index);
s_vec = sim_struct_neq(1).metric_array(:,sharpness_index);
cw_filter = cw_vec >= 5.2 & s_vec >0;

[~,mi] = nanmax(sim_struct_neq(1).metric_array(:,decision_rate_index).*cw_filter);

inputParams = sim_struct_neq(1).rate_array(mi,:);
paramCell = mat2cell(inputParams,size(inputParams,1),ones(1,size(inputParams,2)));
RNum = RSymFun(paramCell{:});

numerical_precision = 10;
state_probs = steadyStateVecFunction(paramCell{:}); 
RNumFlux = RNum.*state_probs;
RNumFlux = RNumFlux / max(RNumFlux(:)) * 10;

RNumGraph = RNum;
% RNumGraph(RNumFlux<=1e-4) = 0;

g1 = digraph(RNumGraph>0);

RNumGraphSmall = RNum; 
RNumGraphSmall (RNumFlux>1e-4) = 0;

g2 = digraph(RNumGraphSmall>0);

close all
figure; 

hold on
% lw = log(RNumFlux(RNum>0)) - min(log(RNumFlux(RNum>0)))+1e-6;
p1 = plot(g1,'LineWidth',RNumFlux(RNum>0));
p1.MarkerSize = 50*state_probs;
p1.XData = [1 2 2 1 0 0];
p1.YData = [1 1 0 0 0 1 ];
%% 
tic
[simInfoNeq, sim_struct_neq] = param_sweep_multi([sharp_spec_index spec_alt_index],sweep_options{:},...
                                          'half_max_flag',true,'wrongFactorConcentration',1,'equilibrium_flag',false);                                        
toc                                                          
%%

[sim_info_spec_eq, sim_struct_spec_eq] = param_sweep_multi([sharpness_index spec_index],sweep_options{:},...
                                          'half_max_flag',true,'wrongFactorConcentration',1,'equilibrium_flag',true);
%%                                        
[sim_info_spec_neq, sim_struct_spec_neq] = param_sweep_multi([sharp_right_norm_index spec_index],sweep_options{:},...
                                          'half_max_flag',true,'wrongFactorConcentration',1,'equilibrium_flag',false,'testingFlag',1); 
                                        
%%                                        
                           
sharp_ind = find(sim_struct_spec_neq(2).metric_array(:,sharpness_index)>=.48&sim_struct_spec_neq(2).metric_array(:,spec_index)>=0,1);      
sharp_rates = sim_struct_spec_neq(2).rate_array(sharp_ind,:);
valMatSharp = [1 sharp_rates];
valCellSharp = mat2cell(valMatSharp,size(valMatSharp,1),ones(1,size(valMatSharp,2)));
state_probs_sharp = steadyStateVecFunction(valCellSharp{:})

spec_ind = find(sim_struct_spec_neq(2).metric_array(:,sharpness_index)>=.22&sim_struct_spec_neq(2).metric_array(:,spec_index)>=1.99,1);                                                                                
spec_rates = sim_struct_spec_neq(2).rate_array(spec_ind,:);                                        
valMatSpec = [1 spec_rates];
valCellSpec = mat2cell(valMatSpec,size(valMatSpec,1),ones(1,size(valMatSpec,2)));                                        
state_probs_spec = steadyStateVecFunction(valCellSpec{:})                                        
                                        
                                        
                                        
%%        
cr_vec = logspace(0,5,10);
info_cell_eq = cell(size(cr_vec));
info_cell_neq = cell(size(cr_vec));

for i = 1:length(cr_vec)
    tic
    
    [sim_info_spec_eq, sim_struct_spec_eq] = param_sweep_multi([affinity_index decision_rate_index],sweep_options{:},...
                                              'half_max_flag',true,'wrongFactorConcentration',cr_vec(i),'equilibrium_flag',true);

    [sim_info_spec_neq, sim_struct_spec_neq] = param_sweep_multi([affinity_index decision_rate_index],sweep_options{:},...
                                          'half_max_flag',true,'wrongFactorConcentration',cr_vec(i),'equilibrium_flag',false); 
                                        
    % find best performing eq and noneq networks
    metric_array_eq = vertcat(sim_struct_spec_eq.metric_array);   
    info_vec_eq = metric_array_eq(:,decision_rate_index);
    [~,info_mi_eq] = nanmax(info_vec_eq);
    info_cell_eq{i} = info_vec_eq(info_mi_eq);
    
    metric_array_neq = vertcat(sim_struct_spec_neq.metric_array);   
    info_vec_neq = metric_array_neq(:,decision_rate_index);
    [~,info_mi_neq] = nanmax(info_vec_neq);
    info_cell_neq{i} = info_vec_neq(info_mi_neq);

    toc
end
%%
id_cell_eq = cell(size(info_cell_neq));
id_cell_neq = cell(size(info_cell_neq));

for i = 1:length(id_cell_eq)
    id_cell_eq{i} = repelem(cr_vec(i),length(info_cell_eq{i}))';
    id_cell_neq{i} = repelem(cr_vec(i),length(info_cell_neq{i}))';
end

close all
figure;
scatter(vertcat(id_cell_eq{:}),vertcat(info_cell_eq{:}))
hold on
scatter(vertcat(id_cell_neq{:}),vertcat(info_cell_neq{:}))
%%
% [sim_info_neq1, sim_struct_neq1] = param_sweep_multi([switch_noise_index spec_index],sweep_options{:},...
%                                                             'wrongFactorConcentration',100,'equilibrium_flag',false);
                                                          
% [sim_info_eq, sim_struct_eq] = param_sweep_multi([sharpness_index spec_index],sweep_options{:},...
%                                                             'half_max_flag',false,'equilibrium_flag',true);                                                          
%                                                           
% [sim_info_eq, sim_struct_eq] = param_sweep_multi([flux_index sharp_spec_index],sweep_options{:},...
%                                                             'half_max_flag',false,'equilibrium_flag',true);                                                          
                                                          
toc                                                          

% [sim_info_neq_half, sim_struct_neq_half] = param_sweep_multi([flux_index decision_rate_index],sweep_options{:},'half_max_flag',true);
% toc
% 
% % set save name
% save_name_flux = ['param_sweep_results_' metric_names{flux_index} '_' ...
%     metric_names{decision_rate_index}];
% save([OutPath save_name_flux '_eq0.mat'],'sim_struct_neq','-v7.3')
% save([OutPath save_name_flux 'info_eq0.mat'],'sim_info_neq')
% 
% save([OutPath save_name_flux '_half_eq0.mat'],'sim_struct_neq','-v7.3')
% save([OutPath save_name_flux '_half_info_eq0.mat'],'sim_info_neq')

%% %%%%%%%%%%%%%%%% info rate vs "affinity" (i.e. log(koff/kon) ) %%%%%%%%%
tic
[simInfoNeq, sim_struct_neq]  = param_sweep_multi([affinity_index decision_rate_index],sweep_options{:},'equilibrium_flag',false);

[sim_info_eq, sim_struct_eq]  = param_sweep_multi([affinity_index decision_rate_index],sweep_options{:},'equilibrium_flag',true);
% sim_struct_neq = param_sweep_multi([affinity_index decision_rate_index],sweep_options{:},'equilibrium_flag',false,'n_seeds',50);
toc

% set save name
% save_name_affinity = ['param_sweep_results_' metric_names{affinity_index} '_' ...
%     metric_names{decision_rate_index}];
% save([OutPath save_name_affinity '_eq0.mat'],'sim_struct_neq','-v7.3')
% save([OutPath save_name_affinity 'info_eq0.mat'],'sim_info_neq')
% save([OutPath save_name_affinity '_eq1.mat'],'sim_struct_eq','-v7.3')
% save([OutPath save_name_affinity 'info_eq1.mat'],'sim_info_eq')

%% %%%%%%%%%%%%%%%% decision time vs "affinity" (i.e. log(koff/kon) ) %%%%%%%%%
tic
[simInfoNeq, sim_struct_neq] = param_sweep_multi([affinity_index decision_time_index],sweep_options{:},'equilibrium_flag',false);

[sim_info_eq, sim_struct_eq] = param_sweep_multi([affinity_index decision_time_index],sweep_options{:},'equilibrium_flag',true);
% sim_struct_neq = param_sweep_multi([affinity_index decision_rate_index],sweep_options{:},'equilibrium_flag',false,'n_seeds',50);
toc

% set save name
% save_name_affinity = ['param_sweep_results_' metric_names{affinity_index} '_' ...
%     metric_names{decision_time_index}];
% save([OutPath save_name_affinity '_eq0.mat'],'sim_struct_neq','-v7.3')
% save([OutPath save_name_affinity 'info_eq0.mat'],'sim_info_neq')
% save([OutPath save_name_affinity '_eq1.mat'],'sim_struct_eq','-v7.3')
% save([OutPath save_name_affinity 'info_eq1.mat'],'sim_info_eq')


%% %%%%%%%%%%%%%%%% sharpness vs precision %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
[simInfoNeq, sim_struct_neq] = param_sweep_multi([sharpness_index precision_index],sweep_options{:},...
                        'equilibrium_flag',false,'half_max_flag',false);
toc

tic
[sim_info_eq, sim_struct_eq] = param_sweep_multi([sharpness_index precision_index],sweep_options{:},...
                        'equilibrium_flag',true,'half_max_flag',false);
toc

% set save name
save_name_tradeoff = ['param_sweep_results_' metric_names{sharpness_index} '_' ...
                                              metric_names{precision_index}];
save([OutPath save_name_tradeoff '_eq0.mat'],'sim_struct_neq','-v7.3')
save([OutPath save_name_tradeoff 'info_eq0.mat'],'sim_info_neq')
save([OutPath save_name_tradeoff '_eq1.mat'],'sim_struct_eq','-v7.3')
save([OutPath save_name_tradeoff 'info_eq1.mat'],'sim_info_eq')

%% %%%%%%%%%%%%%%%% sharpness vs switching noise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_indices = [sharpness_index switch_noise_index];
tic
[simInfoNeq, sim_struct_neq] = param_sweep_multi(plot_indices, sweep_options{:},...
                        'equilibrium_flag',false,'half_max_flag',true);
toc

tic
[sim_info_eq, sim_struct_eq] = param_sweep_multi(plot_indices,sweep_options{:},...
                        'equilibrium_flag',true,'half_max_flag',true);
toc

% set save name
save_name_tradeoff = ['param_sweep_results_' metric_names{plot_indices(1)} '_' ...
                      metric_names{plot_indices(2)}];
save([OutPath save_name_tradeoff '_eq0.mat'],'sim_struct_neq','-v7.3')
save([OutPath save_name_tradeoff 'info_eq0.mat'],'sim_info_neq')
save([OutPath save_name_tradeoff '_eq1.mat'],'sim_struct_eq','-v7.3')
save([OutPath save_name_tradeoff 'info_eq1.mat'],'sim_info_eq')

%% %%%%%%%%%%% sharpness vs precision (with half-max constraint) %%%%%%%%%%
plot_indices = [sharpness_index switch_noise_index];

tic
[simInfoNeq, sim_struct_neq] = param_sweep_multi(plot_indices,sweep_options{:},...
                        'equilibrium_flag',false,'half_max_flag',true);
toc

tic
[sim_info_eq, sim_struct_eq] = param_sweep_multi(plot_indices,sweep_options{:},...
                        'equilibrium_flag',true,'half_max_flag',true);
toc

% set save name
save_name_tradeoff = ['param_sweep_results_' metric_names{plot_indices(1)} '_' ...
                      metric_names{plot_indices(2)}];
save([OutPath save_name_tradeoff '_half_eq0.mat'],'sim_struct_neq','-v7.3')
save([OutPath save_name_tradeoff '_halfinfo_eq0.mat'],'sim_info_neq')
save([OutPath save_name_tradeoff '_half_eq1.mat'],'sim_struct_eq','-v7.3')
save([OutPath save_name_tradeoff '_halfinfo_eq1.mat'],'sim_info_eq')
