% script to simulate sharp and precise circuit dynamics
clear 
close all
addpath(genpath('../utilities/'))

% set basic parameters
nStates = 4;
rate_bounds = repmat([-8 ; 4],1,9); % constrain transition rate magnitude
[~,~,metric_names] = calculateMetricsSym_v2([]);

% make sure we're linked to the appropriate function subfolder
functionPath = ['../utilities/metricFunctions/n' num2str(nStates) '_OR/'];
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));

% define save path
OutPath = ['../../out/bivariate_parameter_sweeps_n' num2str(nStates) filesep];
mkdir(OutPath);
                         

% get index of useful metrics
rate_index = find(strcmp(metric_names,'Production Rate'));
sharpness_norm_index = find(strcmp(metric_names,'SharpnessNormed'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
precision_index = find(strcmp(metric_names,'Precision'));
precision_index_norm = find(strcmp(metric_names,'PrecisionNorm'));
precision_poisson_index = find(strcmp(metric_names,'PrecisionPoisson'));
pp_normed_index = find(strcmp(metric_names,'PrecPoissonNormed'));
ir_index = find(strcmp(metric_names,'DecisionRateNorm'));
tau_index = find(strcmp(metric_names,'CycleTime'));



% set sim options
sweep_options = {'n_seeds',10,'n_iters_max',50,'nStates',nStates,'numCalcFlag',0};

%% %%%%%%%%%%%%%%%% obtain optimal networks %%%%%%%%%%%%%%%%%%%%
TauCycleLimit = 100;
                                      
[simInfoNeq, sim_struct_neq] = param_sweep_multi_v3([sharpness_index precision_index],functionPath,sweep_options{:},...
                                          'half_max_flag',true,'equilibrium_flag',false, 'TauCycleLimit',TauCycleLimit);   
                                        
%% Find optimal sharp an precise networks
[max_sharpness,sharp_i] = nanmax(sim_struct_neq(1).metric_array(:,sharpness_index));
[max_precision,prec_i] = nanmax(exp(sim_struct_neq(1).metric_array(:,precision_index)));
[max_ir,ir_i] = nanmax(sim_struct_neq(1).metric_array(:,ir_index));

sharp_rates = sim_struct_neq(1).rate_array(sharp_i,:);
precise_rates = sim_struct_neq(1).rate_array(prec_i,:);
ir_rates = sim_struct_neq(1).rate_array(ir_i,:);
%%
infoCell = mat2cell(sharp_rates,size(ir_rates,1),ones(1,size(ir_rates,2)));       
R_matrix_info = RSymFun(infoCell{:});

%% Simulate
t_sim = 1e4; % number of time points
dT = 50;
time_grid = 0:dT:t_sim;
initiation_rate = 1/3; % Pol II per second

n_sim = 100; % number of independent simulations to run
n_reps = 1e2; % number of replicates over which to calculate variance
state_options = 1:4;
pd_states = [3,4]; % specify states that make mRNA
init_rate_vec = zeros(size(state_options));
init_rate_vec(pd_states) = initiation_rate;
% basic rate parameters
paramBounds = repmat([-3 ; -1],1,8); % log10 of bounds on transition rates 
c_bounds = [-1 1];


% generate lists of rate and concentration parameters
rng(123)

% calculate steady state stats
sharpCell = mat2cell(sharp_rates,size(sharp_rates,1),ones(1,size(sharp_rates,2)));       
R_matrix_sharp = RSymFun(sharpCell{:});    
[V,D] = eig(R_matrix_sharp);
[~,mi] = max(real(diag(D)));
ss = V(:,mi)/sum(V(:,mi));
jump_weight_array = R_matrix_sharp;
jump_weight_array(eye(4)==1) = 0;

% extract state dwell vector
dwell_time_vec = -diag(R_matrix_sharp)'.^-1;       

% initilize array to store cumulative mRNA
total_mRNA_array = NaN(length(time_grid),n_reps);
total_mRNA_array_poisson = NaN(length(time_grid),n_reps);
on_dwell_time_array = NaN(1,n_reps);
off_dwell_time_array = NaN(1,n_reps);

for rep = 1%:n_reps

  % initialize vector to store results
  jump_time_vec = [];
  off_dwell_times = [];
  on_dwell_times = [];
  initiation_count_vec = [];
  state_vec = [randsample(state_options,1,true,ss)];    
  total_time = 0;   
  act_time = 0;
  while total_time < t_sim     

    % extract state info
    state_curr = state_vec(end);
    active_flag_curr = ismember(state_curr,pd_states);

    % calculate jump time means
    tau = dwell_time_vec(state_curr);    

    % randomly draw next jump time
    dt = exprnd(tau);

    % randomly draw the number of initiation events
    n_pol_II = poissrnd(init_rate_vec(state_curr)*dt);

    total_time = total_time + dt;
    act_time = act_time + dt;

    % randomly draw destination state            
    next_state = randsample(state_options,1,true,jump_weight_array(:,state_curr));      
    active_flag_next = ismember(next_state,pd_states);

    % assign states
    state_vec(end+1) = next_state;
    jump_time_vec(end+1) = dt;    
    initiation_count_vec(end+1) = n_pol_II;

    % update ON/OFF dwell time vecs if appropriate
    if active_flag_next && ~active_flag_curr
        off_dwell_times(end+1) = act_time;
        act_time = 0;
    elseif ~active_flag_next && active_flag_curr
        on_dwell_times(end+1) = act_time;
        act_time = 0;
    end
  end  

  % downsample and save
  active_state_vec = double(ismember(state_vec(1:end-1),pd_states));
  cumulative_time = cumsum(jump_time_vec);
  cumulative_mRNA = cumsum(jump_time_vec.*active_state_vec);
  cumulative_mRNA_poisson = cumsum(initiation_count_vec);
  total_mRNA_array(:,rep) = interp1(cumulative_time,cumulative_mRNA,time_grid);
  total_mRNA_array_poisson(:,rep) = interp1(cumulative_time,cumulative_mRNA_poisson,time_grid);

  % calculate average dwell times in ON and OFF states
  on_dwell_time_array(rep) = mean(on_dwell_times);
  off_dwell_time_array(rep) = mean(off_dwell_times);
end  

% record average values
sim_struct(n).r_mean = mean(total_mRNA_array(end,:))/time_grid(end);
sim_struct(n).r_var = var(total_mRNA_array(end,:))/time_grid(end);
sim_struct(n).r_mean_poisson = mean(total_mRNA_array_poisson(end,:))/time_grid(end);
sim_struct(n).r_var_poisson = var(total_mRNA_array_poisson(end,:))/time_grid(end);
sim_struct(n).mRNA_array = total_mRNA_array;

sim_struct(n).tau_on = mean(on_dwell_time_array);
sim_struct(n).tau_off = mean(off_dwell_time_array);
sim_struct(n).tau_cycle = mean(on_dwell_time_array+off_dwell_time_array);
