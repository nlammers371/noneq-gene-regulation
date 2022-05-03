% This script is intended to run stochastic simulations to verify 
% theoretical expressions for intrinsic noise for the full four state network

clear
close all
addpath(genpath('../../utilities'))


OutPath = '../../out/validation/';
mkdir(OutPath);

%%%%%%%%%%%%%%%%%
% specify core simulation parameters
%%%%%%%%%%%%%%%%%
t_sim = 1e4; % number of time points
dT = 50;
time_grid = 0:dT:t_sim;

n_sim = 100; % number of independent simulations to run
n_reps = 1e2;
state_options = 1:4;
pd_states = [3,4]; % specify states that make mRNA
% basic rate parameters
rate_bounds = [-2 2]; % log10 of bounds on transition rates 
c_bounds = [-1 1];

dT_init = .1; % use constant time scale to generate correlated noise and resample as needed

% generate lists of rate and concentration parameters
rng(123)
rate_array = 10.^(rand(n_sim,8)*(rate_bounds(2)-rate_bounds(1))-rate_bounds(2));
c_vec = 10.^(rand(1,n_sim)*(c_bounds(2)-c_bounds(1))-c_bounds(2));

% initialize parpool
myCluster = parcluster('local');
NumWorkers = myCluster.NumWorkers;
p = gcp('nocreate');
if isempty(p)
  parpool(ceil(NumWorkers/2));
end

% define anonymous function to generate transition rate matrix
R = @(k,r,c) [-k(3)-r(4)  c*k(2)             0               r(1); 
             r(4)      -r(3)-c*k(2)         k(1)              0
              0           r(3)        -c*r(2)-k(1)           k(4)
             k(3)          0              c*r(2)        -r(1)-k(4) ];        
  
% initialize arrays
sim_struct = struct;

% initialize waitbar
h = waitbar(0,'Simulating transcription network dynamics...');
% iterate through time points  
for n = 1:n_sim
  waitbar(n/n_sim,h);  
  
  % extract params
  rate_vec = rate_array(n,:);
  mu_c = c_vec(n);     
  
  % get predictions
  sim_struct(n).r_predicted = fourStateProduction_v4(rate_vec,mu_c);
  sim_struct(n).var_predicted = fourStateVariance_v4(rate_vec,mu_c);
  
  % calculate steady state stats
  R_matrix = R(rate_vec(1:4),rate_vec(5:8),mu_c);    
  [V,D] = eig(R_matrix);
  [~,mi] = max(real(diag(D)));
  ss = V(:,mi)/sum(V(:,mi));
  jump_weight_array = R_matrix;
  jump_weight_array(eye(4)==1) = 0;
  
  % extract state dwell vector
  dwell_time_vec = -diag(R_matrix)'.^-1;       
  
  % initilize array to store cumulative mRNA
  total_mRNA_array = NaN(length(time_grid),n_reps);
  
  parfor rep = 1:n_reps
    
    % initialize vector to store results
    jump_time_vec = [];
    state_vec = [randsample(state_options,1,true,ss)];    
    total_time = 0;   
    tic
    while total_time < t_sim     
      
      % extract state info
      state_curr = state_vec(end);

      % calculate jump time means
      tau = dwell_time_vec(state_curr);    

      % randomly draw jump times
      dt = exprnd(tau);
      total_time = total_time + dt;
      
      % randomly draw destination state            
      next_state = randsample(state_options,1,true,jump_weight_array(:,state_curr));      

      % assigne states
      state_vec(end+1) = next_state;
      jump_time_vec(end+1) = dt;    
    end  
    
    % downsample and save
    cumulative_time = cumsum(jump_time_vec);
    cumulative_mRNA = cumsum(jump_time_vec.*double(ismember(state_vec(1:end-1),pd_states)));
    total_mRNA_array(:,rep) = interp1(cumulative_time,cumulative_mRNA,time_grid);
  end  
  
  % record average values
  sim_struct(n).r_mean = mean(total_mRNA_array(end,:))/time_grid(end);
  sim_struct(n).r_var = var(total_mRNA_array(end,:))/time_grid(end);
  sim_struct(n).mRNA_array = total_mRNA_array;
end

delete(h);

% save results
save([OutPath 'mean_var_sim_validation.mat'],'sim_struct')


