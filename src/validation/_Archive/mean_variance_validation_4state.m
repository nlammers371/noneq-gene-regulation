% script to cofirim that analytic expressions for occupancy and variance of
% 4 state network agree with results of stochastic simulations

clear
close all
addpath('../../utilities')
% define parameters
rate_range = 1;
rng(20); % for consistency
rate_vec = 10.^(rand(1,8)*2*rate_range - rate_range);
% rate_vec = [0.9151    0.6160    1.0000    1.0000    0.0029    0.9324    0.0022    0.3016]*1e4;
c_val = 1;

% define rate matrix 
R =     [-rate_vec(3)-rate_vec(end)       c_val*rate_vec(2)               0                            rate_vec(5); 
             rate_vec(end)            -c_val*rate_vec(2)-rate_vec(7)     rate_vec(1)                           0
              0                              rate_vec(7)     -rate_vec(1)-c_val*rate_vec(6)          rate_vec(4)
             rate_vec(3)                          0                    c_val*rate_vec(6)        -rate_vec(4)-rate_vec(5) ];
         
% calculate predicted production rates and variance
pd_rate_predicted = fourStateProductionGeneral(rate_vec,c_val);
pd_var_predicted = fourStateVarianceGeneral(rate_vec,c_val);

% calculate full SS vec
[V,D] = eig(R);
[~,mi] = max(real(diag(D)));
ss = V(:,mi)/sum(V(:,mi))

%% simulation parameters
T = 500;
n_sim = 100;
t_sim = 0:0.1:T;
% simulate trajectories of independent replicates of stocahstic 
% transcriptional networks

occupancy_vec = NaN(1,n_sim);
production_array = NaN(T,n_sim);
production = [0 0 1 1];
state_options = 1:4;

% iterate
for n = 1:n_sim
    if mod(n,10) == 0
        disp(n)
    end
    
    state_vec = [randsample(state_options,1,true,ss)];
    jump_vec = [];   
    production_vec = [production(state_vec(1))];
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
    occupancy_vec(n) = sum(jump_vec(ismember(state_vec(1:end-1),[3,4])));
    pd_rs = interp1(cumsum([0 jump_vec]),cumsum(production_vec),t_sim);
    production_array(1:numel(pd_rs),n) = pd_rs;
end

%%
mean(occupancy_vec) / T
var(occupancy_vec) / T
cm = jet(128);

figPath = '../../fig/validation/';
mkdir(figPath)
close all

mean_fig = figure;
hold on
mu_exp = nanmean(production_array,2);
plot(t_sim,mu_exp,'LineWidth',1);
plot(t_sim,pd_rate_predicted*t_sim,'--','LineWidth',1);
xlabel('time')
ylabel('mean accumulated mRNA (au)')
grid on
set(gca,'Fontsize',14)
legend('simulation','analytic prediction','Location','northwest')
saveas(mean_fig,[figPath 'mean_analytic_vs_sim.png'])

var_fig = figure;
hold on
var_exp = nanvar(production_array,[],2);
plot(t_sim(1:end-100),var_exp(1:end-100),'LineWidth',1);
plot(t_sim,pd_var_predicted*t_sim,'--')
xlabel('time')
ylabel('variance in accumulated mRNA (au)')
grid on
set(gca,'Fontsize',14)
legend('analytic prediction','simulation','Location','northwest')
saveas(var_fig,[figPath 'var_analytic_vs_sim.png'])



