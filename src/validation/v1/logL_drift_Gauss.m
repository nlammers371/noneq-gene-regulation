% script to validate analytic expression for rate of log likelihood divergence
clear
close all
addpath('../../utilities')

% define parameters
c1 = 1.1;
c0 = 1;
activator_flag = 1;

options = logspace(-2,1);

r_vec = randsample(options,4);
k_vec = randsample(options,4);

R1 = [-c1*k_vec(1)-r_vec(4)       r_vec(1)                  0                    k_vec(4); 
       c1*k_vec(1)           -r_vec(1)-k_vec(2)          r_vec(2)                  0
        0                          k_vec(2)          -r_vec(2)-k_vec(3)        c1*r_vec(3)
       r_vec(4)                       0                    k_vec(3)        -c1*r_vec(3)-k_vec(4) ];
   
R0 = [-c0*k_vec(1)-r_vec(4)       r_vec(1)                  0                    k_vec(4); 
       c0*k_vec(1)           -r_vec(1)-k_vec(2)          r_vec(2)                  0
        0                          k_vec(2)          -r_vec(2)-k_vec(3)        c0*r_vec(3)
       r_vec(4)                       0                    k_vec(3)        -c0*r_vec(3)-k_vec(4) ];  


[V1,D1] = eig(R1);
[~,mi] = max(real(diag(D1)));
ss1 = V1(:,mi)/sum(V1(:,mi));

[V0,D0] = eig(R0);
[~,mi] = max(real(diag(D0)));
ss0 = V0(:,mi)/sum(V0(:,mi));


m1 = fourStateProduction([k_vec r_vec],c1,activator_flag);
m0 = fourStateProduction([k_vec r_vec],c0,activator_flag);

s1 = sqrt(fourStateVariance([k_vec r_vec],c1,activator_flag));
s0 = sqrt(fourStateVariance([k_vec r_vec],c0,activator_flag));
v1 = s1^2;
v0 = s0^2;
%
% calculate expected drift and diffusion rates
V = .25*((m1-m0)^2*(s0^2+s1^2))/(s0^2*s1^2)
D = .25*((m1-m0)^2*(s0^6+s1^6))/(s0^4*s1^4)
% calcualte logL expectation
K = log(100);
B = V*K/D;
tau = K/(2*V*sinh(B)) * (exp(B) + exp(-B) - 2)

% simulation parameters

T = round(tau);
T = max(1000,round(tau));
T = min(10000,T);
n_sim = 50;
t_grid = 1:T;

% calculate t values for K values
log_vec = linspace(0,K,1000);

b_vec = V*log_vec./D;
tau_vec = log_vec./(2*V*sinh(b_vec)) .* (exp(b_vec) + exp(-b_vec) - 2);


% simulate trajectories and track movement of log likelihood ratio over
% time
occupancy_vec = NaN(1,n_sim);
production_array = NaN(T,n_sim);
logL_array = NaN(T,n_sim);
production = [0 0 1 0];
state_options = 1:4;
% iterate
for n = 1:n_sim
    if mod(n,10) == 0
        disp(n)
    end    
    jump_vec = [];       
    production_vec = [0];
    total_time = 0;    
    if mod(n,2) == 0
        R = R0;
        state_vec = [randsample(state_options,1,true,ss0)];
    else
        R = R1;
        state_vec = [randsample(state_options,1,true,ss1)];
    end
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
    pd_rs = interp1(cumsum([0 jump_vec]),cumsum(production_vec),1:T);
    production_array(:,n) = pd_rs;
    % calculate log likelihood ratio
    t_vec = 1:T;
    logL_array(:,n) = abs(log(sqrt(v0/v1)) +  .5*(m0.*t_vec-pd_rs).^2./v0./t_vec...
        -.5*(m1.*t_vec-pd_rs).^2./v1./t_vec);
end

rate_fig = figure;
hold on
plot(nanmean(production_array,2));
plot((.5*m1+.5*m0)*(1:T));

var_fig = figure;
hold on
index_vec = 1:n_sim;
even_vec = 2:2:n_sim;
plot(.5*nanstd(production_array(:,even_vec),[],2)+.5*nanstd(production_array(:,~ismember(index_vec,even_vec)),[],2));
plot(sqrt((.5*v1+.5*v0)*(1:T)));

logL_fig = figure;
hold on
logL_vec = nanmean(logL_array,2);
% delta = logL_vec(2300) - V*2300;
plot(1:T,logL_vec);
plot(tau_vec,log_vec);
plot(V*t_vec+logL_vec(1))
xlim([0 T])

