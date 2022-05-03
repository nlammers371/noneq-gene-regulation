% script to validate analytic expression for variance of general 4 state
% Markov chain using stochastic simulations
clear
close all
addpath('../utilities')
% define parameters
c = 3;
activator_flag = 1;

r1 = rand()*10;
r2 = rand()*10;
r3 = rand()*10;
r4 = rand()*10;

k1 = rand()*10;
k2 = rand()*10;
k3 = rand()*10;
k4 = rand()*10;

r_vec = [r1,r2,r3,r4];
k_vec = [k1,k2,k3,k4];

R = [-c*k1-r4       r1             0          k4; 
       c*k1       -r1-k2          r2          0
        0         k2          -r2-k3        c*r3
       r4          0            k3        -c*r3-k4 ];
tic       
[V,D] = eig(R);
toc
[~,mi] = max(real(diag(D)));
ss = V(:,mi)/sum(V(:,mi));
rateAlt = ss(3)
ratePredicted = ...
    fourStateProduction([k_vec r_vec],c,activator_flag)

varPredicted = ...
    fourStateVariance([k_vec r_vec],c,activator_flag)


% simulation parameters
T = 1000;
n_sim = 50;
t_grid = 1:T;
% simulate trajectories and track movement of log likelihood ratio over
% time
occupancy_vec = NaN(1,n_sim);
production_array = NaN(T,n_sim);
logL_cell = cell(1,n_sim);
production = [0 0 1 0];
state_options = 1:4;
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
        next_state = randsample(options,1,true,R(options,current_state));
        total_time = total_time + dt;
        if total_time < T
            jump_vec = [jump_vec dt];
            state_vec = [state_vec next_state];
            production_vec = [production_vec dt*production(current_state)];                        
        end
    end
    occupancy_vec(n) = sum(jump_vec(state_vec(1:end-1)==3));
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

% snr_fig = figure;
% hold on
% plot(snr_vec_sim)
% plot(snr_vec_pd)
% legend('simulation','analytic solution','Location','northwest')
% grid on
