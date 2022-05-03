% script to simulate sharp and precise circuit dynamics
clear 
close all
addpath(genpath('../utilities/'))

% set basic parameters
nStates = 4;
paramBounds = repmat([-3 ; 3],1,9); % constrain transition rate magnitude
[~,~,metric_names] = calculateMetricsSym_v2([]);

% make sure we're linked to the appropriate function subfolder
functionPath = ['../utilities/metricFunctions/n' num2str(nStates) '_OR/'];
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));
                         
% define save path
OutPath = ['../../out/illustrative_bursting_simulations' filesep];
mkdir(OutPath);

% get index of useful metrics
rate_index = find(strcmp(metric_names,'Production Rate'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
precision_index = find(strcmp(metric_names,'Precision'));
ir_index = find(strcmp(metric_names,'DecisionRateNorm'));


% set sim options
sweep_options = {'n_seeds',10,'n_iters_max',50,'nStates',nStates,'numCalcFlag',0};

% %%%%%%%%%%%%%%%% obtain optimal networks %%%%%%%%%%%%%%%%%%%%
TauCycleTime = 1;
                                      
[sim_info_neq, sim_struct_neq] = param_sweep_multi_v3([sharpness_index precision_index],functionPath,sweep_options{:},...
                                                    'equilibrium_flag',false, 'TauCycleTime',TauCycleTime,'paramBounds',paramBounds);  
[sim_info_eq, sim_struct_eq] = param_sweep_multi_v3([sharpness_index precision_index],functionPath,sweep_options{:},...
                                                    'equilibrium_flag',true, 'TauCycleTime',TauCycleTime,'paramBounds',paramBounds);                                                    
                                        
% Find optimal sharp an precise networks
rate_bounds = [0.495 0.505];
rate_filter_neq = sim_struct_neq.metric_array(:,rate_index)>rate_bounds(1) & sim_struct_neq.metric_array(:,rate_index)<rate_bounds(2) & sim_struct_neq.metric_array(:,sharpness_index) > 0;
rate_filter_eq = sim_struct_eq.metric_array(:,rate_index)>rate_bounds(1) & sim_struct_eq.metric_array(:,rate_index)<rate_bounds(2) & sim_struct_eq.metric_array(:,sharpness_index) > 0;

[max_ir_neq,ir_i_neq] = nanmax(sim_struct_neq.metric_array(:,ir_index).*(1*rate_filter_neq));
[max_ir_eq,ir_i_eq] = nanmax(sim_struct_eq.metric_array(:,ir_index).*(1*rate_filter_eq));

ir_rates_neq = sim_struct_neq.rate_array(ir_i_neq,:);
ir_rates_eq = sim_struct_eq.rate_array(ir_i_eq,:);


% Simulate
t_sim = 2e2; % number of time burst cycles to simulate
dT = 1e-2;
time_grid = 0:dT:t_sim;
initiation_rate = 1; % unit-less production rate

n_reps = 1e2; % number of replicates over which to calculate variance
state_options = 1:4;
pd_states = [3,4]; % specify states that make mRNA

% initialize structure to track simulation results
c0 = 0.83; % note that we're making hte difference larger for illustrative purposes
c1 = 1.17;
sim_struct = struct;
sim_struct(1).name = 'nonequilibrium (c0)';
sim_struct(1).rate_vec = ir_rates_neq;
sim_struct(1).rate_vec(1) = c0;
sim_struct(1).c = c0;
sim_struct(2).name = 'nonequilibrium (c1)';
sim_struct(2).rate_vec = ir_rates_neq;
sim_struct(2).rate_vec(1) = c1;
sim_struct(2).c = c1;
sim_struct(3).name = 'equilibrium (c0)';
sim_struct(3).rate_vec = ir_rates_eq;
sim_struct(3).rate_vec(1) = c0;
sim_struct(3).c = c0;
sim_struct(4).name = 'equilibrium (c1)';
sim_struct(4).rate_vec = ir_rates_eq;
sim_struct(4).rate_vec(1) = c1;
sim_struct(4).c = c1;

% initialize random number generator
rng(123)

for s = 1:length(sim_struct)
  
    sim_rates = sim_struct(s).rate_vec;
    
%     wb = waitbar(0,['Simulating gene circuits for system type: ' sim_struct(s).name]);
    
    % calculate key simulation stats for neq system
    info_cell = mat2cell(sim_rates,size(sim_rates,1),ones(1,size(sim_rates,2)));       
    R_matrix = RSymFun(info_cell{:});
    [V,D] = eig(R_matrix);
    [~,mi] = max(real(diag(D)));
    ss_vec = V(:,mi)/sum(V(:,mi));
    jump_weight_array = R_matrix;
    jump_weight_array(eye(4)==1) = 0;   

    % extract state dwell vector
    dwell_time_vec = -diag(R_matrix)'.^-1;               

    % initilize array to store cumulative mRNA
    total_mRNA_array = NaN(length(time_grid),n_reps);
    promoter_state_array = NaN(length(time_grid),n_reps);

    for rep = 1:n_reps
%       waitbar(rep/n_reps,wb);
      
      % initialize vector to store results
      jump_time_vec = [];       
      state_vec = [randsample(state_options,1,true,ss_vec)];      
      total_time = 0;   

      while total_time < t_sim     

        % extract state info
        state_curr = state_vec(end);
        active_flag_curr = ismember(state_curr,pd_states);

        % calculate jump time means
        tau = dwell_time_vec(state_curr);    

        % randomly draw next jump time
        dt = exprnd(tau);
        total_time = total_time + dt;

        % randomly draw destination state            
        next_state = randsample(state_options,1,true,jump_weight_array(:,state_curr));      

        % assign states
        state_vec(end+1) = next_state;
        jump_time_vec(end+1) = dt;    

      end  

      % downsample and save
      active_state_vec = double(ismember(state_vec(1:end),pd_states));
      cumulative_time = [0 cumsum(jump_time_vec)];
      cumulative_mRNA = cumsum(jump_time_vec.*active_state_vec(1:end-1));  
      total_mRNA_array(:,rep) = interp1([cumulative_time],[0 cumulative_mRNA],time_grid);  
      promoter_state_array(:,rep) = interp1([cumulative_time],[active_state_vec],time_grid,'previous');  

    end      
    sim_struct(s).total_mRNA_array = total_mRNA_array;
    sim_struct(s).promoter_state_array = promoter_state_array;
    sim_struct(s).time_grid = time_grid;
%     delete(wb);
end

save([OutPath 'sim_struct.mat'], 'sim_struct')
