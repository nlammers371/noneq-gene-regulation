clear 
close all
addpath('../utilities')
outPath = '../../out/param_optimization/';
mkdir(outPath);
figPath = '../../fig/param_optimization/';
mkdir(figPath);

% generate input profile data
x_vec = linspace(0,1,25);
% f_diff = 1.2;
scale_factor = log(10);%log(f_diff) / (x_vec(2)-x_vec(1));
c_prefactor = exp(scale_factor/2);
c_profile = c_prefactor*exp(-x_vec.*scale_factor);

% iterate through neighboring points and calculate expected decision times
% under different constraints

% params
edge_opt_cell = {4,5,6,7,8,5:8};
activator_flag = 0;
K = log(100);
lb = 1e-4;
ub = 1e1;
n_reps = 10;
% initialize arrays
decision_time_array = NaN(numel(c_profile)-1,n_reps,numel(edge_opt_cell));
opt_rate_array = NaN(numel(c_profile)-1,8,n_reps,numel(edge_opt_cell));
options = optimoptions('fmincon','Display','off');
for r = 1:numel(edge_opt_cell)
    opt_indices = unique([1:4 edge_opt_cell{r}]);
    nr = numel(opt_indices);
    tic
    c1 = c_profile(1); % this remains constant
    parfor c = 1:numel(c_profile)-1       
        c0 = c_profile(c+1);
        for n = 1:n_reps
            % set optimization routine params
            options = logspace(-2,1);
            rate_vec_init = randsample(options,nr);
            % define caller function
            rate_fun_call = @(rate_vec) rate_opt_fun(rate_vec,opt_indices,K,c1,c0,activator_flag);
            fit = fmincon(rate_fun_call,rate_vec_init,[],[],[],[],repelem(lb,nr),repelem(ub,nr));
            % record
            [tau_opt,rate_vec_opt] = rate_opt_fun(fit,opt_indices,K,c1,c0,activator_flag);
            opt_rate_array(c,:,n,r) = rate_vec_opt;
            decision_time_array(c,n,r) = tau_opt;
        end
    end
    toc
end

sim_struct = struct;
sim_struct.edge_opt_vec = edge_opt_cell;
sim_struct.c_profile = c_profile;
sim_struct.rate_bounds = [lb ub];
im_struct.activator_flag = activator_flag;
sim_struct.logL_rate = K;
sim_struct.decision_time_array = decision_time_array;
sim_struct.opt_rate_array = opt_rate_array;
% save 
save([outPath 'tau_fold_dependence_act' num2str(activator_flag) '.mat'],'sim_struct');