% script to explore decision convergence ideas via simulation. Assumes 4
% state system in which two transition rates are concentration-dependent
clear
close all
% define parameters
phi = .1;
koff = 10;
L0 = 1;
L1 = 10;
L = L1;
kp = .1;
km = .1;
beta = 100;

k_vec1 = [L1*phi,kp,koff,beta*km];
r_vec1 = [koff,km,L1*phi,kp];

k_vec0 = [L0*phi,kp,koff,beta*km];
r_vec0 = [koff,km,L0*phi,kp];

k_vec = @(x)[x*phi,kp,koff,beta*km];
r_vec = @(x)[koff,km,x*phi,kp];

R = @(k,r) [-k(1)-r(4)   r(1)             0               k(4); 
             k(1)        -r(1)-k(2)      r(2)              0
              0           k(2)         -r(2)-k(3)         r(3)
             r(4)          0               k(3)        -r(3)-k(4) ];
         
R_actual = R(k_vec(L),r_vec(L));    
[V,D] = eig(R_actual);
[~,mi] = max(diag(D));
ss = V(:,mi)/sum(V(:,mi));

tau_actual = 1 / (-R_actual(3,3) * ss(3));
tr = tau_actual - 1 / -R_actual(3,3);
k_actual = 1/tr;

R1 = R(k_vec(L1),r_vec(L1));
[V1,D1] = eig(R1);
[~,mi] = max(diag(D1));
ss1 = V1(:,mi)/sum(V1(:,mi));
tr1 = -R1(3,3)^-1 * (1-ss1(3)) / ss1(3);
k1 = 1/tr1;

R0 = R(k_vec(L0),r_vec(L0));
[V0,D0] = eig(R0);
[~,mi] = max(diag(D0));
ss0 = V0(:,mi)/sum(V0(:,mi));
tr0 = -R0(3,3)^-1 * (1-ss0(3)) / ss0(3);
k0 = 1/tr0;
 
% simulation parameters
T = 15000;
n_sim = 1;
% simulate trajectories and track movement of log likelihood ratio over
% time
state_cell = cell(1,n_sim);
jump_time_cell = cell(1,n_sim);
jump_time_approx_cell = cell(1,n_sim);
log_ratio_full_cell = cell(1,n_sim);
log_ratio_approx_cell = cell(1,n_sim);
state_options = 1:4;
% iterate
for n = 1:n_sim
    state_vec = [3];
    jump_vec = [];
    jump_vec_approx = [];
    log_ratio_full_vec = [0];    
    log_ratio_approx_vec = [0];    
    total_time = 0;    
    while total_time < T
        current_state = state_vec(end);
        % draw jump time
        dt = exprnd(1/-R_actual(current_state,current_state),1);
        options = state_options(state_options~=current_state);
        next_state = randsample(options,1,true,R_actual(options,current_state));
        total_time = total_time + dt;
        if total_time < T
            jump_vec = [jump_vec dt];
            state_vec = [state_vec next_state];
            % full likelihood
            if (next_state == 2 && current_state == 1) || (next_state == 3 && current_state == 4)
                l_new = log_ratio_full_vec(end) + log(L1/L0)-dt*(L1-L0)*phi; 
                log_ratio_full_vec = [log_ratio_full_vec l_new];
            else
                log_ratio_full_vec = [log_ratio_full_vec log_ratio_full_vec(end)];
            end
            % approximate likelihood
            if next_state == 3 
                indices = find(state_vec==3,2,'last');
                dt_eff = sum(jump_vec(indices(end-1)+1:end));
                jump_vec_approx = [jump_vec_approx dt_eff];
                l_new = log_ratio_approx_vec(end) + log(k1/k0)-dt_eff*(k1-k0); 
                log_ratio_approx_vec = [log_ratio_approx_vec l_new];
            else
                log_ratio_approx_vec = [log_ratio_approx_vec log_ratio_approx_vec(end)];
            end
        end
    end
    log_ratio_full_cell{n} = log_ratio_full_vec;
    log_ratio_approx_cell{n} = log_ratio_approx_vec;
    state_cell{n} = state_vec;
    jump_time_cell{n} = jump_vec;
    jump_time_approx_cell{n} = jump_vec_approx;
end


%%% compare simulation to analytical prediction
sim_log_dt = (log_ratio_full_vec(end) - log_ratio_full_vec(1)) / T
sim_log_approx_dt = (log_ratio_approx_vec(end) - log_ratio_approx_vec(1)) / T
koff_eff = koff+km;
analytic_log_dt = (k_actual*koff_eff) / (k_actual+koff_eff) * ((k0 - k1)/k_actual +  log(k1/k0))


%%% make figure
t_axis = cumsum([jump_time_cell{1}]);
sim_comparison = figure;
hold on
plot(t_axis,log_ratio_approx_cell{1}(2:end))
plot(t_axis,t_axis*analytic_log_dt)
legend('simulation','prediction','Location','northwest')

tau_sim = nanmean(jump_vec_approx)
tau_calc = 1/k_actual


