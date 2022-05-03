% script to explore decision convergence ideas via simulation. Assumes 4
% state system in which two transition rates are concentration-dependent
clear
close all
% define par
ameters
phi = 1;
koff = 10;
L0 = 1;
L1 = 1;
L = L1;
kp = 1;
km = 1;
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
 
% simulation parameters
T = 5000;
n_sim = 100;
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
        dt = exprnd(1/-R_actual(current_state,current_state),1);
        options = state_options(state_options~=current_state);
        next_state = randsample(options,1,true,R_actual(options,current_state));
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

%%
mean(occupancy_vec) / T
var(occupancy_vec) / T
cm = jet(128);

mean_fig = figure;
mu = nanmean(production_array(1:T-1,:),2);
plot(production_array(1:T-10,:),'LineWidth',1);
xlabel('time')
ylabel('accumulated mRNA (au)')
grid on
figPath = '../fig/simulation_results/';
mkdir(figPath)
saveas(mean_fig,[figPath 'mean_mt.png'])

sigma = std(production_array(1:T-1),[],2);
ub = mu + sigma;
lb = mu - sigma;

spread_fig = figure;
hold on
fill([1:T-1 fliplr(1:T-1)],[ub' fliplr(lb')],cm(30,:),'FaceAlpha',.2,'EdgeAlpha',0)
fill([1:T-1 fliplr(1:T-1)],[ub' fliplr(lb')]/2,cm(120,:),'FaceAlpha',.2,'EdgeAlpha',0)

p1 = plot(1:T-1,mu,'Color',cm(30,:));
p2 = plot(1:T-1,mu/2,'Color',cm(120,:));

xlabel('time')
ylabel('accumulated mRNA (au)')
grid on
legend([p1,p2],'region 1','region 0','Location','northwest')
ylim([0 max(ub)])
saveas(spread_fig,[figPath 'spread_mt.png'])

