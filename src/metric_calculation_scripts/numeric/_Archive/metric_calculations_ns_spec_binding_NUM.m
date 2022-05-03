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

n_bs_vec = 1:5;
n_g_vec = 1:5;
n_state_vec = (n_bs_vec+1).*(n_g_vec+1)';

for m = 1:length(n_g_vec)
    for n = 1:length(n_bs_vec)    
        nStates = n_state_vec(m,n);

        savePath = [writePath 'n' sprintf('%03d',nStates) '_s' sprintf('%02d',n_bs_vec(n)) '_ns00_g' sprintf('%02d',n_g_vec(m))  filesep];
        mkdir(savePath);
        addpath(savePath)

        % call script to extend model
        [networkInfo, RSym, RSymFull] = makeActivatorGeneralNetwork(networkInfoInit,n_bs_vec(n),n_g_vec(m));

        save([savePath 'networkInfo'],'networkInfo');     
    %     clear networkInfo
        matlabFunction(RSym,'File',[savePath 'RSymFun'],'Optimize',true,...
                                    'Vars',networkInfo.sweepVarList);             
        
        if nStates <= 16
            matlabFunction(RSymFull,'File',[savePath 'RSymFunFull'],'Optimize',true,...
                                  'Vars',networkInfo.sweepVarList);
        end

    end        
end