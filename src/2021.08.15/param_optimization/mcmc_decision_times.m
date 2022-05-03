clear 
close all
addpath('../utilities')
outPath = '../../out/mcmc_optimization/';
mkdir(outPath);
figPath = '../../fig/mcmc_optimization/';
mkdir(figPath);

% generate input profile data
x_vec = linspace(0,1,25);
f_diff = 1.1;
scale_factor = log(f_diff) / (x_vec(2)-x_vec(1));
c_prefactor = exp(scale_factor/2);
c_profile = c_prefactor*exp(-x_vec.*scale_factor);

% iterate through neighboring points and calculate expected decision times
% under different constraints

% params
rate_key = [1 2 3 4; 
            7 8 5 6]; 
edge_opt_vec = 4:8;
activator_flag = 0;
K = log(100);
rate_bounds = [1e-4,1e1];
prop_sigma = .2;
n_samples = 500;
ref_indices = 1:8;
temp = 100;
burn_in = 200;
options = optimoptions('fmincon','Display','off');
% initialize arrays
decision_time_array = NaN(numel(c_profile)-1,n_samples,numel(edge_opt_vec));
opt_rate_array = NaN(numel(c_profile)-1,8,n_samples,numel(edge_opt_vec));
for r = 2%1:numel(edge_opt_vec)
    opt_indices = unique([1:4 edge_opt_vec(r)]);
    nr = numel(opt_indices);
    stable_indices = find(~ismember(rate_key(2,:),opt_indices));
    for c = 1:numel(c_profile)-1
        c1 = c_profile(c);
        c0 = c_profile(c+1);
        % set optimization routine params
        init_tau_vec = NaN(1,10);
        init_rate_array = NaN(10,8);
        % define caller function
        rate_fun_call = @(rate_vec) rate_opt_fun(rate_vec,opt_indices,K,c1,c0,activator_flag);
%         for i = 1:10
%             options = logspace(-1,1);
%             prop_rate_vec = randsample(options,nr);                
%             fit = fmincon(rate_fun_call,prop_rate_vec,[],[],[],[],repelem(rate_bounds(1),nr),...
%                 repelem(rate_bounds(2),nr));
%             [init_tau_vec(i),init_rate_array(i,:)] = rate_opt_fun(fit,opt_indices,K,c1,c0,activator_flag);
%         end
%         [current_tau, mi] = min(init_tau_vec);
%         current_rate_vec = init_rate_array(mi,:);       
        % perform MCMC sampling        
        for n = 1:n_samples            
            % record current values
            opt_rate_array(c,:,n,r) = current_rate_vec;
            decision_time_array(c,n,r) = current_tau;
            % sample new rates
            new_rates = sample_rates_general(current_rate_vec(opt_indices),rate_bounds,prop_sigma);  
            new_rate_vec = NaN(1,8);
            new_rate_vec(opt_indices) = new_rates;
            new_rate_vec(rate_key(2,stable_indices)) = new_rates(rate_key(1,stable_indices));                      
            % calculate new decision time
            new_tau = rate_fun_call(new_rates);
            % MH move
            score = exp((current_tau-new_tau)/temp);
            update = score > rand();
            if update
                current_rate_vec = new_rate_vec;
                current_tau = new_tau;
            end            
        end              
    end
end

sim_struct = struct;
sim_struct.edge_opt_vec = edge_opt_vec;
sim_struct.c_profile = c_profile;
sim_struct.rate_bounds = rate_bounds;
im_struct.activator_flag = activator_flag;
sim_struct.logL_rate = K;
sim_struct.decision_time_array = decision_time_array;
sim_struct.opt_rate_array = opt_rate_array;
% save 
save([outPath 'rate_optimization_results_act' num2str(activator_flag) '.mat'],'sim_struct');

%% make figure
decison_time_fig = figure;
hold on
cmap = brewermap(128,'Spectral');
colormap(cmap);
for i = 1:numel(edge_opt_vec)
    if i == 1
        ss = 's';
    else
        ss = 'o';
    end
    scatter(x_vec(1:end-1)+.05,nanmean(decision_time_array(:,burn_in:end,i),2),30,ss,'MarkerFaceColor',cmap(20*(i-1)+5,:)...
        ,'MarkerEdgeAlpha',0);
end
grid on
set(gca,'yscale','log')
ylabel('decision time')
xlabel('position')
