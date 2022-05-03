% script to explore decision convergence ideas via simulation. Assumes 2
% state system in which *both* transition rates are concentration-dependent
clear
close all
% define parameters
phi = 2;
theta = 1;
L0 = 1;
L1 = 2;
L = L1;
R = [-phi*L theta*L; phi*L -theta*L];
% simulation parameters
T = 10000;
n_sim = 1;
states = [1 2];
% simulate trajectories and track movement of log likelihood ratio over
% time
state_cell = cell(1,n_sim);
jump_time_cell = cell(1,n_sim);
log_ratio_cell = cell(1,n_sim);

% iterate
for n = 1:n_sim
    state_vec = [1];
    jump_vec = [];
    log_ratio_vec = [0];    
    total_time = 0;
    
    while total_time < T
        current_state = state_vec(end);
        % draw jump time
        dt = exprnd(1/-R(current_state,current_state),1);
        total_time = total_time + dt;
        if total_time < T
            jump_vec = [jump_vec dt];
            state_vec = [state_vec states(states~=current_state)];
   
            l_new = log_ratio_vec(end) + log(L1/L0)-dt*(L1-L0)*R(state_vec(end),current_state)/L; 
            log_ratio_vec = [log_ratio_vec l_new];

        end
    end
    log_ratio_cell{n} = log_ratio_vec;
    state_cell{n} = state_vec;
    jump_time_cell{n} = jump_vec;
end


%%% compare simulation to analytical prediction
sim_log_dt = (log_ratio_vec(end) - log_ratio_vec(1)) / T

analytic_log_dt = (theta*phi) / (theta+phi) * 2 * (L0 - L1 +  L * log(L1/L0))


%%% make figure
t_axis = cumsum([0 jump_time_cell{1}]);
sim_comparison = figure;
hold on
plot(t_axis,log_ratio_cell{1})
plot(t_axis,t_axis*analytic_log_dt)
legend('simulation','prediction','Location','northwest')
