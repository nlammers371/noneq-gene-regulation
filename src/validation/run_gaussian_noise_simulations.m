% This script is intended to run stochastic simulations to verify 
% theoretical expressions for intrinsic noise for the full four state network
clear
close all
addpath(genpath('../utilities'))

DropboxFolder = 'S:\Nick\Dropbox\Nonequilibrium\Nick\SweepOutput';
OutPath = [DropboxFolder filesep 'appendices_v4' filesep];
mkdir(OutPath);

% set path to approapriate functions
% functionPath = '../utilities/metricFunctions/n4_OR';
functionPath = '../utilities/metricFunctions/symbolic/n004_s01_ns00_g01';
rmpath(genpath('../utilities/metricFunctions/'));
addpath(genpath(functionPath));

%%%%%%%%%%%%%%%%%
% specify core simulation parameters
%%%%%%%%%%%%%%%%%
t_sim = 5e3; % duration of simulation (in burst cycles)
dT = 5;
time_grid = logspace(-2,log10(t_sim),1000);%0:dT:t_sim;
% initiation_rate = 1/3; % Pol II per second

n_sim = 500; % number of independent simulations to run
n_traces = 1e2; % number of replicates over which to calculate variance
n_init_factor = 100;
state_options = 1:4;
pd_states = [3,4]; % specify states that make mRNA
a_bounds = [0.02 0.98];
init_rate_vec = zeros(size(state_options));
init_rate_vec(pd_states) = 1; % work in units of r0
% basic rate parameters
paramBounds = repmat([-2 ; 2],1,8);
paramBounds(1,5) = 0; % ensures activation
paramBounds(2,6) = 0; % ensures activation
c_bounds = [-1 1];

% generate lists of rate and concentration parameters
rng(123)

% generate arrays of upper and lower bounds
lb_array = repmat(paramBounds(1,:),n_sim*n_init_factor,1);
ub_array = repmat(paramBounds(2,:),n_sim*n_init_factor,1);

% draw sample rates
rate_array = reshape(10.^(trandn(lb_array,ub_array)),[],size(paramBounds,2)); % function to sample from truncated normal distribution
c_vec = 10.^(rand(1,n_sim*n_init_factor)*(c_bounds(2)-c_bounds(1))-c_bounds(2)); % activator values

% define anonymous function to generate transition rate matrix
param_array = [c_vec' rate_array];
paramCell = mat2cell(param_array,size(param_array,1),ones(1,size(param_array,2)));

% identify systems that fall within prescribed bounds
ProductionRateArray = productionRateFunction(paramCell{:});

% calculate predicted metrics
TauOnArray = TauONFunction(paramCell{:});
TauOffArray = TauOFFFunction(paramCell{:});
TauCycleArray = TauOffArray+TauOnArray;  

% convert to units of cycle time
param_array(:,2:5) = param_array(:,2:5) .* TauCycleArray;
paramCell = mat2cell(param_array,size(param_array,1),ones(1,size(param_array,2)));
VarianceArray = intrinsicVarianceFunction(paramCell{:});

option_flags = ProductionRateArray>=a_bounds(1) & ProductionRateArray<=a_bounds(2)...
    &all(log10(param_array(:,2:end))>=paramBounds(1,:),2)&all(log10(param_array(:,2:end))<=paramBounds(2,:),2);

if sum(option_flags) < n_sim
    error('Not enough systems meet production rate requirement. Try generating more initial options.')
else    
    sim_indices = randsample(find(option_flags),n_sim,'false');    
    ProductionRateArray = ProductionRateArray(sim_indices,:);
%     rate_array = rate_array(sim_indices,:);
%     c_vec = c_vec(sim_indices);
    param_array = param_array(sim_indices,:);%[c_vec' rate_array];
    paramCell = mat2cell(param_array,size(param_array,1),ones(1,size(param_array,2)));
end

% initialize parpool
myCluster = parcluster('local');
NumWorkers = myCluster.NumWorkers;
p = gcp('nocreate');
if isempty(p)
  parpool(28);%ceil(NumWorkers/1.5));
end      
  
% initialize arrays
Gaussian_noise_struct = struct;
% 
% %% initialize waitbar
% Gaussian_noise_struct = struct;
% rate_array = rate_array(empty_inds,:);
% mu_c = c_vec(empty_inds);
% ProductionRateArray = ProductionRateArray(empty_inds);
% VarianceArray = VarianceArray(empty_inds);
% param_array = param_array(empty_inds,:);
% %%
% h = waitbar(0,'Simulating transcription network dynamics...');
% iterate through time points  
parfor n = 1:length(ProductionRateArray)
%   n = empty_inds(n)
  tic
%   waitbar(n/n_sim,h);  
  
  % extract params
%   rate_vec = rate_array(n,:);
  mu_c = c_vec(n);     
  
  % get predictions
  Gaussian_noise_struct(n).r_predicted = ProductionRateArray(n);
  Gaussian_noise_struct(n).var_predicted = VarianceArray(n);
%   sim_struct(n).r_predicted_poisson = sim_struct(n).r_predicted*initiation_rate;
%   sim_struct(n).var_predicted_poisson = sim_struct(n).var_predicted*initiation_rate^2 +  sim_struct(n).r_predicted_poisson;
  
  %   sim_struct(n).tau_on_predicted = TauOnArray(n);
%   sim_struct(n).tau_off_predicted = TauOffArray(n);
%   sim_struct(n).tau_cycle_predicted = TauCycleArray(n);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Calculate quantities needed for simulation  
  tempCell = mat2cell(param_array(n,:),size(param_array(n,:),1),ones(1,size(param_array(n,:),2)));       
  R_matrix = RSymFun(tempCell{:});    
  [V,D] = eig(R_matrix);
  [~,mi] = max(real(diag(D)));
  ss = V(:,mi)/sum(V(:,mi));
  jump_weight_array = R_matrix;
  jump_weight_array(eye(4)==1) = 0;
  
  % extract state dwell vector
  dwell_time_vec = -diag(R_matrix)'.^-1;       
  
  % initilize array to store cumulative mRNA
  total_mRNA_array = NaN(length(time_grid),n_traces);
  
  for rep = 1:n_traces
    
    % initialize vector to store results
    jump_time_vec = [];
%     off_dwell_times = [];
%     on_dwell_times = [];
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
      dt = exprnd(tau) + 1e-30; % prevents rare numerical errors
      
      % randomly draw the number of initiation events
%       n_pol_II = poissrnd(init_rate_vec(state_curr)*dt);
      
      total_time = total_time + dt;
      act_time = act_time + dt;
      
      % randomly draw destination state            
      next_state = randsample(state_options,1,true,jump_weight_array(:,state_curr));      
      active_flag_next = ismember(next_state,pd_states);
      
      % assign states
      state_vec(end+1) = next_state;
      jump_time_vec(end+1) = dt;    
      initiation_count_vec(end+1) = init_rate_vec(state_curr)*dt;
      
    end  
    
    % downsample and save
    active_state_vec = double(ismember(state_vec(1:end-1),pd_states));
    try
        cumulative_time = [realmin cumsum(jump_time_vec)];
    %     cumulative_mRNA = [0 cumsum(jump_time_vec.*active_state_vec)];
        cumulative_mRNA = [0 cumsum(initiation_count_vec)];
        total_mRNA_array(:,rep) = interp1(cumulative_time,cumulative_mRNA,time_grid);
    catch
        warning('Error occurred during interpolation')
    end
%     total_mRNA_array_poisson(:,rep) = interp1(cumulative_time,cumulative_mRNA_poisson,time_grid);
    
    % calculate average dwell times in ON and OFF states
%     on_dwell_time_array(rep) = mean(on_dwell_times);
%     off_dwell_time_array(rep) = mean(off_dwell_times);
  end  
  
  % record average values
  Gaussian_noise_struct(n).r_mean = nanmean(total_mRNA_array(end,:))/time_grid(end);
  Gaussian_noise_struct(n).r_var = nanvar(total_mRNA_array(end,:))/time_grid(end);
  Gaussian_noise_struct(n).r_mean_vec = nanmean(total_mRNA_array,2)./time_grid';
  Gaussian_noise_struct(n).r_var_vec = nanvar(total_mRNA_array,[],2)./time_grid';
%   sim_struct(n).r_mean_poisson = mean(total_mRNA_array_poisson(end,:))/time_grid(end);
%   sim_struct(n).r_var_poisson = var(total_mRNA_array_poisson(end,:))/time_grid(end);
  Gaussian_noise_struct(n).mRNA_array = total_mRNA_array;
  Gaussian_noise_struct(n).time_grid = time_grid;
  
%   sim_struct(n).tau_on = mean(on_dwell_time_array);
%   sim_struct(n).tau_off = mean(off_dwell_time_array);
%   sim_struct(n).tau_cycle = mean(on_dwell_time_array+off_dwell_time_array);
  
  % perform normality test
  total_mRNA_array_norm = (total_mRNA_array - nanmean(total_mRNA_array,2)) ./ nanstd(total_mRNA_array,[],2);
%   Gaussian_noise_struct(n).total_mRNA_array_norm = total_mRNA_array_norm;
  Gaussian_noise_struct(n).p_vec = NaN(size(time_grid));
  Gaussian_noise_struct(n).gauss_flag = NaN(size(time_grid));
  for t = 1:length(time_grid)
      try
          [Gaussian_noise_struct(n).gauss_flag(t),Gaussian_noise_struct(n).p_vec(t)] = kstest(total_mRNA_array_norm(t,:));
      catch
          % do nothing
      end
  end
      
  toc 
end

% delete(h);
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
% %%
% fnames = fieldnames(Gaussian_noise_struct);
% 
% for i = 1:length(empty_inds)
%     for f = 1:length(fnames)
%         Gaussian_noise_struct_orig(empty_inds(i)).(fnames{f}) = Gaussian_noise_struct(i).(fnames{f});
%     end
% end

save([OutPath 'Gaussian_noise_struct.mat'],'Gaussian_noise_struct', '-v7.3')
%%
Gaussian_noise_struct = rmfield(Gaussian_noise_struct,'total_mRNA_array_norm');
%%
for i = 1:length(Gaussian_noise_struct)
    total_mRNA_array = Gaussian_noise_struct(i).mRNA_array;
    total_mRNA_array_ds = interp1(time_grid,total_mRNA_array,0:5:t_sim);
    Gaussian_noise_struct(i).mRNA_array = total_mRNA_array_ds;
    Gaussian_noise_struct(i).time_grid = 0:5:t_sim;
end    

