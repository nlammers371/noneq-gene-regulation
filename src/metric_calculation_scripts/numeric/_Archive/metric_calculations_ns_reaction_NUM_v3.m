% This generates increasingly complex architecture files for use in
% numerical parameter sweeps. We assume 1 non-specific reaction
% and an increasing number of specific activator binding reactions, which
% are assumed to be identical for simplicity
clear
close all

% generate path to save metric functions 
subfolderName = 'numeric';
writePath = handlePathOptions(subfolderName);

% first, generate baseline 6 state network that serves as kernel for all
% architectures
networkInfoInit = generateBaselineSixStateNetwork;

% generate generic prefix for function directories 
folder_string = 's01_ns00_g';
%

%% generate systems with 2-10 binding sites
% as we add additional reaction, this will multiply the state space by
% a factor of 2, since each state from smaller "base" network can coincide
% with the new factor being (1) engaged or (2) disengaged

n_ns_vec = 1:5;
n_state_vec = 2*(n_ns_vec+1);


for n = 2:length(n_ns_vec)    
    nStates = n_state_vec(n);
    
    savePath = [writePath 'n' sprintf('%03d',nStates) '_' folder_string sprintf('%02d',n_ns_vec(n)) filesep];
    mkdir(savePath);
    addpath(savePath)
    
    % call script to extend model
    [networkInfo, RSym, RSymFull] = makeGenericNSNetwork(networkInfoInit,n_ns_vec(n));  
    save([savePath 'networkInfo'],'networkInfo');     
    
%     clear networkInfo
    matlabFunction(RSym,'File',[savePath 'RSymFun'],'Optimize',true,...
                                'Vars',networkInfo.sweepVarList);             
    
    matlabFunction(RSymFull,'File',[savePath 'RSymFunFull'],'Optimize',true,...
                          'Vars',networkInfo.sweepVarList);
    
end        
