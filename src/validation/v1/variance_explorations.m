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
N = 1000; % number of visits (not time)
n_sim = 100;
t_grid = 1:N;
% simulate trajectories and track movement of log likelihood ratio over
% time
occupancy_vec = NaN(1,n_sim);
production_vec = NaN(1,n_sim);
logL_cell = cell(1,n_sim);
production = [0 0 0 1];
mRNA_state = find(production==1);
state_options = 1:4;
% iterate
for n = 1:n_sim
    if mod(n,10) == 0
        disp(n)
    end
    state_vec = [1];         
    total_mRNA = 0;
    total_time = 0;
    total_visits = 0;    
    while total_visits < N
        current_state = state_vec(end);
        % draw jump time
        dt = exprnd(1/-R(current_state,current_state),1);
        options = state_options(state_options~=current_state);
        next_state = randsample(options,1,true,R(options,current_state));
        if current_state == mRNA_state
            total_visits = total_visits + 1;
            total_mRNA = total_mRNA + dt;
        end
        total_time = total_time + dt;
        if total_visits < N            
            state_vec = [state_vec next_state];            
        end
    end
    occupancy_vec(n) = total_mRNA / total_time;    
    production_vec(n) = total_mRNA;
end
%%
var_pd = R(mRNA_state,mRNA_state)^-2 * N
var_sim = var(production_vec) 

