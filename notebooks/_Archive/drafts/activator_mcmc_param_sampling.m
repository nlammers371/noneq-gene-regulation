% Script to conduct parameter search for fidelity of forward-driven
% transcriptional network
clear
close all
addpath('./utilities/')
% list of function options that embody different assumptions
function_list = {'FidelityMaxNorm2','FidelityMaxNorm1'};
parameter_list = {[{'kon'},{'km'},{'kp'},{'c1'}], [{'kon'},{'km'},{'kp'},{'c1'},{'koff'}]};
lb_vec = [0 0 0 1+1e-6 0];
% specify simulation parameters
n_iterations = 1e3;
init_vectors = logspace(-3,3);
n_runs = 10;
temperature = .1;
% specify function to use
function_id = 1;
n_params = numel(parameter_list{function_id});
% initialize arrays to store sampling results
fidelity_array = NaN(n_iterations,n_runs);
param_array = NaN(n_iterations,n_params,n_runs);
% conduct sampling
for sim = 1:n_runs
    init = randsample(init_vectors,n_params,true);    
    init(4) = randsample(logspace(0,3),1) + .1;
    param_array(1,:,sim) = init; 
    jump_vec = param_array(1,:,sim) / 10;
    current_objective = eval([function_list{function_id} '(init(1),init(2),init(3),init(4))']);    
    fidelity_array(1,sim) = current_objective;
    for iter = 2:n_iterations
        current_params = param_array(iter-1,:,sim);
        for i = 1:n_params
            % generate jump distribution
            sigma = jump_vec(i);
            cv = current_params(i);
            pd = makedist('Normal','mu',cv,'sigma',sigma');
            t = truncate(pd,lb_vec(i),Inf);
            proposal = random(t,1);                        
            % draw sample
            thresh = rand();
            param_options = [cv proposal];        
            new_params = current_params;
            new_params(i) = proposal;     
            
            new_objective = eval([function_list{function_id} ...
                '(new_params(1),new_params(2),new_params(3),new_params(4))']);
            
            [~, mi] = max([thresh exp(-(current_objective-new_objective)/temperature)]);
            objective_options = [current_objective new_objective];
            current_objective = objective_options(mi);
            current_params(i) = param_options(mi);               
        end
        fidelity_array(iter,sim) = current_objective;
        param_array(iter,:,sim) = current_params;
        disp(['completed ' num2str(iter) ' of ' num2str(n_iterations) ' (' num2str(sim) ' of ' num2str(n_runs) ')']) 
    end
end            
    

%% Simple Scatters to examine relationship between fidelity and various parameters
fidelity_vec = fidelity_array(:);
for i = 1:n_params
    param_name = parameter_list{function_id}{i};
    param_vals = reshape(param_array(:,i,:),1,[]);
    scatter(fidelity_vec,param_vals);
    error('asfa')
end