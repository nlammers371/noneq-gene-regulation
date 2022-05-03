% This generates increasingly complex architecture files for use in
% numerical parameter sweeps. We assume 1-5 non-specific reaction
% and 1-5 specific activator binding reactions all reactions of each type
% are assumed identical for simplicity

clear
close all

% generate path to save metric functions 
subfolderName = 'numeric';
writePath = handlePathOptions(subfolderName);

% first, generate baseline 6 state network that serves as kernel for all
% architectures
networkInfoInit = generateBaselineSixStateNetwork;

% generate generic prefix for function directories 
% folder_string = 'ns00_g01';       
%% generate systems with 2-10 binding sites
% as we add additional reaction, this will multiply the state space by
% a factor of 2, since each state from smaller "base" network can coincide
% with the new factor being (1) engaged or (2) disengaged

n_bs_vec = 1:3;
n_g_vec = 1:3;
n_ns_vec = 1:3;
n_bs_cw = [6 12 20 30 42]/2;
n_state_vec = (n_bs_vec+1).*(n_g_vec+1)';
n_state_vec_cw = (n_bs_cw).*(n_g_vec+1)';
ns_flag = 1;
for i = 2%1:length(n_ns_vec)%0:1
    for m = 2%:length(n_g_vec)
        for n = 1:2%1:length(n_bs_vec)                                                          

            % call script to extend model           
            [networkInfo, RSym, RSymFull] = makeActivatorGeneralNetworkNSWrong(networkInfoInit,n_bs_vec(n),n_ns_vec(i),n_g_vec(m));            

            nStates = size(RSym,1);
            savePath = [writePath 'n' sprintf('%03d',nStates) '_s' sprintf('%02d',n_bs_vec(n)) '_ns' sprintf('%02d',n_ns_vec(i)) ...
                        '_g' sprintf('%02d',n_g_vec(m)) '_cw' num2str(ns_flag) filesep];
                      
            mkdir(savePath);
            addpath(savePath)
            
            save([savePath 'networkInfo'],'networkInfo');             
            matlabFunction(RSym,'File',[savePath 'RSymFun'],'Optimize',true,...
                                        'Vars',networkInfo.sweepVarList);             

            if nStates <= 16
                matlabFunction(RSymFull,'File',[savePath 'RSymFunFull'],'Optimize',false,...
                                      'Vars',networkInfo.sweepVarList);
            end

        end        
    end
end    