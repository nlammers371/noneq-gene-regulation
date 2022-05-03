% This script is intended to verify that theoretical expectations for
% intrinsic and extrinsic noise contributions agree with results from
% full stochastic simulations

clear
close all

FigPath = '../../fig/validation_v2/';
mkdir(FigPath);

% specify core simulation parameters
kon = .001; % seconds
koff = 1; % seconds
mu_c = 10; % average activator concentration
t_sim = 1e5; % # time points
n_sim = 3e2; % number of independent simulations to run
tau_c = 50; % timescale of extrinsic fluctuations
dT = min([.1 .1/(mu_c*kon) .1/koff]); % resolution of simulation (must be a lot smaller than rates)
time_vec = linspace(0,t_sim*dT,t_sim);

sigma_c = 1; % variation in concentration
tau_i = 1/(mu_c*kon+koff); % characteristic time scale

state_options = [1 2];
rate_vec = [kon koff];

% kernel to introduce temporal correlation
corr_kernel = exp(-(cumsum(ones(1,4*tau_c/dT)*dT)-dT)/tau_c);

% initialize array to store results
state_array = NaN(t_sim,n_sim);
state_array(1,:) = 1;


% draw uncorrelated random gaussian array to be used for all
% simulations
%tic
tic
gauss_array = normrnd(0,sigma_c,t_sim+4*numel(corr_kernel),n_sim);
% generate correlated noise vec 
corr_gauss_array_raw = conv2(gauss_array,corr_kernel','full');

% normalize
norm_vec = conv(ones(size(gauss_array(:,1))),corr_kernel,'full');
norm_mat = repmat(norm_vec,1,n_sim);
corr_gauss_array_raw = corr_gauss_array_raw(2*numel(corr_kernel):end-3*numel(corr_kernel),:)./...
    norm_mat(2*numel(corr_kernel):end-3*numel(corr_kernel),:);


% scale up noise to desired level    
signs = sign(corr_gauss_array_raw);
corr_gauss_diffs_new = sqrt((corr_gauss_array_raw.^2).*sigma_c^2./var(corr_gauss_array_raw)).*signs;
corr_gauss_array = corr_gauss_diffs_new + mu_c;
corr_gauss_array(corr_gauss_array<0) = 0;
toc

% conduct simulation
tic
state_array(1,:) = randsample(state_options,n_sim,true); 
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
toc    
%% plot results against prediction
close all
mu_out = mu_c*kon/(mu_c*kon+koff);
tau_cycle = (mu_c*kon+koff) / (mu_c*kon*koff);
sigma_int = 2*kon*mu_c*koff/(kon*mu_c+koff)^3;

sigma_ext =  2*(sigma_c)^2*(koff*kon/(kon*mu_c+koff)^2)^2 * (tau_c + (mu_c*kon+koff) / (mu_c*kon*koff)) ;

sigma_total = sigma_int + sigma_ext;

cum_array = cumsum(state_array-1,1).*dT;
mean_vec = nanmean(cum_array,2);
var_vec = nanvar(cum_array,[],2);


fig = figure;
hold on
plot(time_vec,var_vec);
plot(time_vec,sigma_int*time_vec);
plot(time_vec,sigma_total*time_vec);
xlabel('time (seconds)')
ylabel('noise')
set(gca,'Fontsize',14)
legend('total noise','intrinsic noise','Location','northwest')
% saveas(fig,[FigPath 'example_tau' num2str(tau_cycle) '_tc' num2str(tau_c) '.png'])

