% This script is intended to verify that theoretical expectations for
% intrinsic and extrinsic noise contributions agree with results from
% full stochastic simulations

clear
close all
addpath('../../utilities')


OutPath = '../../../out/validation_v2/';
mkdir(OutPath);
sim_label = '4D_sim_v1';

%%%%%%%%%%%%%%%%%
% specify core simulation parameters
%%%%%%%%%%%%%%%%%
vec_dim = 5;
t_sim = 1e5; % number of time points
ds_factor = 1e3;
n_sim = 5e2; % number of independent simulations to run

% basic rate parameters
mu_c = 10; % average activator concentration
var_c_vec = linspace(0.05*mu_c,0.2*mu_c,vec_dim); % variation in concentration
koff_vec = linspace(1e-3,1,vec_dim); % seconds
kon_vec = koff_vec/mu_c; % seconds

% set timescale of fluctuations
tau_c_vec = linspace(1e0,100,vec_dim); % timescale of extrinsic fluctuations

dT_init = .1; % use constant time scale to generate correlated noise and resample as needed

% generate lists of parameters
kon_sim_list = NaN(1,vec_dim^4);
koff_sim_list = NaN(1,vec_dim^4);
tau_c_sim_list = NaN(1,vec_dim^4);
sigma_c_sim_list = NaN(1,vec_dim^4);
iter = 1;

for sc = 1:vec_dim
    sigma_c = var_c_vec(sc);
    for tc = 1:vec_dim
        tau_c = tau_c_vec(tc);
        for k1 = 1:vec_dim
            kon = kon_vec(k1);
            for k2 = 1:vec_dim
                koff = koff_vec(k2);
                % update lists
                tau_c_sim_list(iter) = tau_c;
                sigma_c_sim_list(iter) = sigma_c;
                koff_sim_list(iter) = koff;
                kon_sim_list(iter) = kon;
                iter = iter + 1;
            end
        end
    end
end

% initialize arrays
system_sim_array = NaN(t_sim/ds_factor,n_sim,iter-1);
c_sim_array = NaN(t_sim/ds_factor,n_sim,iter-1);
mean_r_vec = NaN(1,n_sim);
var_r_vec = NaN(1,n_sim);
var_c_vec = NaN(1,n_sim);

% initialize waitbar 
WaitMessage = parfor_wait(length(kon_sim_list), 'Waitbar', true);

try
    parpool(12);
catch
    warning('parpool already exists')
end

% iterate through list of simulation conditions
parfor n = 1:length(kon_sim_list)
  % update waitbar    
  WaitMessage.Send;
  % extract params
  kon = kon_sim_list(n);
  koff = koff_sim_list(n);
  sigma_c = sigma_c_sim_list(n);
  tau_c = tau_c_sim_list(n);
 
  dT = round(t_sim*min([dT_init 0.1*tau_c 0.1/(mu_c*kon) 0.1/koff]))/t_sim; % resolution of simulation (must be a lot smaller than rates)
  time_vec = linspace(0,t_sim*dT,t_sim)';
  time_vec_ds = linspace(0,t_sim*dT,t_sim/ds_factor)';
  
  % calculate mean transcription rate 
  r_mean = mu_c*kon/(mu_c*kon+koff);
  
  % kernel to introduce temporal correlation
  corr_kernel = exp(-(cumsum(ones(1,4*tau_c/dT_init)*dT_init)-dT_init)/tau_c);

  % draw uncorrelated random gaussian array to be used for all
  % simulations  
  gauss_array = normrnd(0,sigma_c,t_sim*dT/dT_init+4*numel(corr_kernel),n_sim);

  % generate correlated noise vec 
  corr_gauss_array_raw = conv2(gauss_array,corr_kernel','full');
  
  % normalize
  norm_vec = conv(ones(size(gauss_array(:,1))),corr_kernel,'full');
  norm_mat = repmat(norm_vec,1,n_sim);
  corr_gauss_array_raw = corr_gauss_array_raw(2*numel(corr_kernel):end-3*numel(corr_kernel),:)./...
      norm_mat(2*numel(corr_kernel):end-3*numel(corr_kernel),:);
  
  % rescale to time tres needed for current simulation
  corr_gauss_array_interp = interp1((1:size(corr_gauss_array_raw,1))*dT_init,corr_gauss_array_raw,(1:t_sim)*dT);


  % scale up noise to desired level    
  signs = sign(corr_gauss_array_interp);
  corr_gauss_diffs_new = sqrt((corr_gauss_array_interp.^2).*sigma_c^2./var(corr_gauss_array_interp)).*signs;
  corr_gauss_array = corr_gauss_diffs_new + mu_c;
  corr_gauss_array(corr_gauss_array<0) = 0; % remove negative values
  
  % initialize array to store results
  state_array = NaN(t_sim,n_sim);
  state_array(1,:) = randsample(1:2,n_sim,true,[1-r_mean r_mean]); 
  for t = 2:t_sim
      % extract state info
      state_curr_vec = state_array(t-1,:);
      swap_vec = double(state_curr_vec==1) + 1;

      % calculate jump time means
      tau_vec = repelem(1/koff,n_sim);
      tau_vec(state_curr_vec==1) = 1./(corr_gauss_array(t,state_curr_vec==1)*kon);

      % randomly draw jump times
      jump_times = exprnd(tau_vec);

      % assigne states
      state_array(t,:) = state_curr_vec;
      state_array(t,jump_times<dT) = swap_vec(jump_times<dT);    
  end  
  % calculate burn-in time
  burn_in = 5*ceil(1/(mu_c*kon+koff) / dT);
  % calculate accumulated mRNA
  cum_array = cumsum(state_array-1,1).*dT;
  % record average values
  mean_r_vec(n) = mean(mean(state_array(burn_in:end,:)-1,2));
  var_r_vec(n) = mean(var(cum_array(burn_in:end,:),[],2)./time_vec(burn_in:end));
  % downsample array 
  cum_array_ds = interp1(time_vec,cum_array,time_vec_ds);
  corr_gauss_array_ds = interp1(time_vec,corr_gauss_array,time_vec_ds);
  % store
  system_sim_array(:,:,n) = cum_array_ds;
  c_sim_array(:,:,n) = corr_gauss_array_ds;
end
WaitMessage.Destroy;
% structure to save simulation output
sim_struct = struct;
% basic characteristics
sim_struct.sim_label = sim_label;
sim_struct.mu_c = mu_c;
sim_struct.vec_dim = vec_dim;
sim_struct.t_sim = t_sim; 
sim_struct.ds_factor = ds_factor;
sim_struct.n_sim = n_sim;
% parameter vectors
sim_struct.tau_c_sim_list = tau_c_sim_list;
sim_struct.sigma_c_sim_list = sigma_c_sim_list;
sim_struct.kon_sim_list = kon_sim_list;
sim_struct.koff_sim_list = koff_sim_list;
% sim results
sim_struct.system_sim_array = system_sim_array;
sim_struct.c_sim_array = c_sim_array;
sim_struct.mean_r_vec = mean_r_vec;
sim_struct.var_r_vec = mean_r_vec;
% save results
save([OutPath sim_label '_simulation_results.mat'],'sim_struct')
