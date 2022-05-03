% script to validate analytic expression for variance of general 4 state
% Markov chain using stochastic simulations
clear
close all
addpath('../utilities')
% define parameters
c = 3;
activator_flag = 1;

r1 = rand()*10;

k1 = rand()*10;

r_vec = [r1];
k_vec = [k1];

R = [-c*k1       r1; 
      c*k1      -r1];
tic       
[V,D] = eig(R);
toc
[~,mi] = max(real(diag(D)));
ss = V(:,mi)/sum(V(:,mi));
ratePredicted = ss(2)
varPredicted = 2*c*k1*r1/(r1+c*k1)^3    

% simulation parameters
T = 500;
n_sim = 100;
t_grid = 1:T;
% simulate trajectories and track movement of log likelihood ratio over
% time
occupancy_vec = NaN(1,n_sim);
production_array = NaN(T,n_sim);
production = [0 1];
state_options = 1:2;
% iterate
for n = 1:n_sim
    if mod(n,10) == 0
        disp(n)
    end
    state_vec = [1];
    jump_vec = [];   
    logL_vec = [0];
    production_vec = [0];
    total_time = 0;    
    while total_time < T
        current_state = state_vec(end);
        % draw jump time
        dt = exprnd(1/-R(current_state,current_state),1);
        options = state_options(state_options~=current_state);
        next_state = options;%randsample(options,1,true,R(options,current_state));
        total_time = total_time + dt;
        if total_time < T
            jump_vec = [jump_vec dt];
            state_vec = [state_vec next_state];
            production_vec = [production_vec dt*production(current_state)];                        
        end
    end
    occupancy_vec(n) = sum(jump_vec(state_vec(1:end-1)==2));
    pd_rs = resample(cumsum(production_vec),cumsum([0 jump_vec]),1);
    production_array(1:numel(pd_rs),n) = pd_rs;
end

gamma_sim = var(occupancy_vec) / T

t_vec = 1:T;
var_vec_pd = varPredicted*t_vec;
production_vec_pd = ratePredicted * t_vec; 
snr_vec_pd = production_vec_pd ./ sqrt(var_vec_pd);

cumulative_occ_array = cumsum(production_array);
var_vec_sim = nanvar(production_array,[],2);
production_vec_sim = nanmean(production_array,2);
snr_vec_sim = production_vec_sim ./ sqrt(var_vec_sim);

cm = jet(128);

mean_rate_fig = figure;
hold on
plot(production_vec_pd)
plot(production_vec_sim)
legend('prediction','simulation')

sigma_rate_fig = figure;
hold on
plot(sqrt(var_vec_pd))
plot(sqrt(var_vec_sim))
legend('prediction','simulation')
