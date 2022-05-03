% script to explore decision convergence ideas via simulation
clear
close all
% define parameters
phi = 1;
v = 1;
L0 = 1;
L1 = 1.01;
L = L1;
R = [-phi*L v; phi*L -v];
% simulation parameters
T = 1000;
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
            if current_state == 2
                log_ratio_vec = [log_ratio_vec log_ratio_vec(end)];
            else                
                l_new = log_ratio_vec(end) + log(L1/L0)-dt*(L1-L0)*phi; 
                log_ratio_vec = [log_ratio_vec l_new];
            end
        end
    end
    log_ratio_cell{n} = log_ratio_vec;
    state_cell{n} = state_vec;
    jump_time_cell{n} = jump_vec;
end


%% compare simulation to analytical prediction
sim_log_dt = (log_ratio_vec(end) - log_ratio_vec(1)) / T

analytic_log_dt = (v*phi) / (v+phi*L) * (L0 - L1 + L * log(L1/L0))
