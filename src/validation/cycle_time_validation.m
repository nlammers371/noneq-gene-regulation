% This script is intended to run stochastic simulations to verify 
% theoretical expressions for burst cycle times for models with two
% activator binding sites or two locus activation steps
clear
close all
addpath(genpath('../utilities'))

% DropboxFolder = 'S:\Nick\Dropbox\Nonequilibrium\Nick\SweepOutput';
% OutPath = [DropboxFolder filesep 'appendices' filesep];
% mkdir(OutPath);

% set basic parameters
n_bs_vec = 1:5;
n_g_vec = 1:5;
ns_flag = 0;

% this contains paths used to address correct functions
functionPathCell = cell(length(n_g_vec),length(n_bs_vec));
saveNameCell = cell(length(n_g_vec),length(n_bs_vec));

% get metric names
[~,~,metric_names] = calculateMetricsNumeric_v3([]);

% generate path to save metric functions 
subfolderName = 'numeric';
sourcePath = handlePathOptions(subfolderName);


for m = 1:length(n_g_vec)
    for n = 1:length(n_bs_vec) 

        % add path to functions        
        readName = ['s' sprintf('%02d',n_bs_vec(n)) '_ns00_g'...
                  sprintf('%02d',n_g_vec(m)) '_cw' num2str(ns_flag)];
        functionPath = [sourcePath readName filesep];
        functionPathCell{m,n} = functionPath;
        saveNameCell{m,n} = readName;

    end
end    

%%%%%%%%%%%%%%%%%
%% specify core simulation parameters
%%%%%%%%%%%%%%%%%
na = 2;
nb = 1;

functionPath = functionPathCell{na,nb};

% make sure we're linked to the appropriate function subfolder
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));

% get metric names
[~,~,metric_names] = calculateMetricsNumeric_v3([]);
rate_index = find(strcmp(metric_names,'ProductionRate'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
tau_b_index = find(strcmp(metric_names,'TauCycle'));

% generate sweep info structure by calling sweep alg with 0 iterations
[sweep_info, ~, ~] = ...
                        param_sweep_multi_v3([rate_index sharpness_index],...
                        functionPath,'n_sim',1,'n_iters_max',0,'useParpool',0);

%%                      
t_sim = 1e4; % duration of simulation (in burst cycles)
dT = 5;
time_grid = 0:dT:t_sim;

n_sim = 1e1; % number of independent simulations to run
n_traces = 1e1; % number of replicates over which to calculate variance
% n_init_factor = 100;
state_options = 1:length(sweep_info.activeStateFilter);
pd_states = find(sweep_info.activeStateFilter); % specify states that make mRNA
% a_bounds = [0.02 0.98];
init_rate_vec = zeros(size(state_options));
init_rate_vec(pd_states) = 1; % work in units of r0

% basic rate parameters
paramBounds = repmat([-2 ; 2],1,length(sweep_info.sweepVarList)-1);
c_bounds = [-1 1];

% generate lists of rate and concentration parameters
rng(123)

% generate arrays of upper and lower bounds
lb_array = repmat(paramBounds(1,:),1,1);
ub_array = repmat(paramBounds(2,:),1,1);

% draw sample rates


% define anonymous function to generate transition rate matrix
% param_array = [c_vec rate_array];
% paramCell = mat2cell(param_array,size(param_array,1),ones(1,size(param_array,2)));

% get predicted cycle times
numerical_precision = 10;
activeStateFilter = sweep_info.activeStateFilter;
sim_struct = struct;
iter = 1;

while iter <= n_sim

    rate_vec = reshape(10.^(trandn(lb_array,ub_array)),[],size(paramBounds,2)); % function to sample from truncated normal distribution
%     c_vec = ones(n_sim,1);
    valCellInit = mat2cell([1 rate_vec],size([1 rate_vec],1),ones(1,size([1 rate_vec],2)));    

    % get rate arrays
    Q_num_cs = RSymFun(valCellInit{:});            
    cs_probs_temp = calculate_ss_num(Q_num_cs,numerical_precision);    

    % calculate cycle time (will use for normalization)
%     used_flag = 1;
    try
        [TauON,TauOFF,TauB] = ...
                calculate_tau_num(Q_num_cs,cs_probs_temp,activeStateFilter,numerical_precision);
        if ~isnan(TauB)            
            sim_struct(iter).TauON = TauON;
            sim_struct(iter).TauOFF = TauOFF;
            sim_struct(iter).TauB = TauB;
            sim_struct(iter).RSym = Q_num_cs;
            sim_struct(iter).ss = cs_probs_temp;
            iter = iter + 1;
        end
    catch
        % do nothing
    end
end            

%%

% initialize parpool
% myCluster = parcluster('local');
% NumWorkers = myCluster.NumWorkers;
% p = gcp('nocreate');
% if isempty(p)
%   parpool(24);%ceil(NumWorkers/1.5));
% end      
  
% calculate predicted metrics
% VarianceArray = intrinsicVarianceFunction(paramCell{:});
% TauOnArray = TauONFunction(paramCell{:});
% TauOffArray = TauOFFFunction(paramCell{:});
% TauCycleArray = TauOffArray+TauOnArray;  
% 
% % convert to units of cycle time
% rate_array(:,2:end) = rate_array(:,2:end) .* TauCycleArray;

% initialize arrays


% initialize waitbar
h = waitbar(0,'Simulating transcription network dynamics...');
% iterate through time points  
for n = 1:n_sim
    tic
    waitbar(n/n_sim,h);  

    % extract params
    rate_vec = rate_array(n,:);
    mu_c = c_vec(n);     

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate quantities needed for simulation  
    tempCell = mat2cell(param_array(n,:),size(param_array(n,:),1),ones(1,size(param_array(n,:),2)));       
    R_matrix = sim_struct(n).RSym;  
    ss = sim_struct(n).ss;
    jump_weight_array = R_matrix;
    jump_weight_array(eye(size(R_matrix,1))==1) = 0;

    % extract state dwell vector
    dwell_time_vec = -diag(R_matrix)'.^-1;       

    off_dwell_time_array = NaN(1,n_traces);
    on_dwell_time_array = NaN(1,n_traces);
    
    % initilize array to store cumulative mRNA
    total_mRNA_array = NaN(length(time_grid),n_traces);

    for rep = 1:n_traces

        % initialize vector to store results
        jump_time_vec = [];
        initiation_count_vec = [];
        off_dwell_times = [];
        on_dwell_times = [];        
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
    %       n_pol_II = poissrnd(init_rate_vec(state_curr)*dt);

          total_time = total_time + dt;
          act_time = act_time + dt;

          % randomly draw destination state            
          next_state = randsample(state_options,1,true,jump_weight_array(:,state_curr));      
          active_flag_next = ismember(next_state,pd_states);

          if active_flag_next==0 && active_flag_curr==1
              on_dwell_times(end+1) = act_time;
              act_time = 0;
          elseif active_flag_next==1 && active_flag_curr==0
              off_dwell_times(end+1) = act_time;
              act_time = 0;
          end
          % assign states
          state_vec(end+1) = next_state;
          jump_time_vec(end+1) = dt;    
          initiation_count_vec(end+1) = init_rate_vec(state_curr)*dt;

        end  

        % downsample and save
        active_state_vec = double(ismember(state_vec(1:end-1),pd_states));
        cumulative_time = [-1e-6 cumsum(jump_time_vec)];    
        cumulative_mRNA = [0 cumsum(initiation_count_vec)];
        total_mRNA_array(:,rep) = interp1(cumulative_time,cumulative_mRNA,time_grid);    

        % calculate average dwell times in ON and OFF states
        on_dwell_time_array(rep) = mean(on_dwell_times);
        off_dwell_time_array(rep) = mean(off_dwell_times);
    end  

    % record average values
%     sim_struct(n).r_mean = mean(total_mRNA_array(end,:))/time_grid(end);
%     sim_struct(n).r_var = var(total_mRNA_array(end,:))/time_grid(end);
%     sim_struct(n).r_mean_vec = mean(total_mRNA_array,2)./time_grid';
%     sim_struct(n).r_var_vec = var(total_mRNA_array,[],2)./time_grid';
  %   sim_struct(n).r_mean_poisson = mean(total_mRNA_array_poisson(end,:))/time_grid(end);
  %   sim_struct(n).r_var_poisson = var(total_mRNA_array_poisson(end,:))/time_grid(end);
%     sim_struct(n).mRNA_array = total_mRNA_array;
%     sim_struct(n).time_grid = time_grid;

    sim_struct(n).tau_on = mean(on_dwell_time_array);
    sim_struct(n).tau_off = mean(off_dwell_time_array);
    sim_struct(n).tau_cycle = mean(on_dwell_time_array+off_dwell_time_array);

%     % perform normality test
%     total_mRNA_array_norm = (total_mRNA_array - mean(total_mRNA_array,2)) ./ std(total_mRNA_array,[],2);
%   %   Gaussian_noise_struct(n).total_mRNA_array_norm = total_mRNA_array_norm;
%     sim_struct(n).p_vec = NaN(size(time_grid));
%     sim_struct(n).gauss_flag = NaN(size(time_grid));
%     for t = 1:length(time_grid)
%         try
%             [sim_struct(n).gauss_flag(t),sim_struct(n).p_vec(t)] = kstest(total_mRNA_array_norm(t,:));
%         catch
%             % do nothing
%         end
%     end

    toc 
end

delete(h);

%%
tau_pd = [sim_struct.TauB];
tau_sim = [sim_struct.tau_cycle];
close all
figure;
scatter(tau_pd,tau_sim)
% %%
% close all
% figure;
% plot(time_grid/sim_struct(n).tau_cycle,sim_struct(n).r_mean_vec*time_grid(end)./time_grid')
% 
% figure;
% hold on
% var_array = [];
% for n = 1:length(sim_struct)
%     var_array(:,n) = sim_struct(n).r_var_vec./ sim_struct(n).var_predicted;
%     plot(time_grid/sim_struct(n).tau_cycle,sim_struct(n).r_var_vec ./ sim_struct(n).var_predicted)
% end
% 
% plot(nanmean(var_array,2),'-k','LineWidth',4)
% save results
% save([OutPath 'Gaussian_noise_struct.mat'],'Gaussian_noise_struct', '-v7.3')
