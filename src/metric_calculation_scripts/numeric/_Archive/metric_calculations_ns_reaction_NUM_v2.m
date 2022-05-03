% This generates increasingly complex architecture files for use in
% numerical parameter sweeps. We assume 1 specific binding site and an
% increasing number of distinct, non-specific reactions 
clear
close all

% generate path to save metric functions 
subfolderName = 'numeric\basic';
writePath = handlePathOptions(subfolderName);

% first, generate baseline 6 state network that serves as kernel for all
% architectures
networkInfoInit = generateBaselineSixStateNetwork;

% generate generic prefix for function directories 
folder_prefix = 's01_ns00_g';

%% generate full 8 network info
% This corresponds to network with one specific activator site and
% two general transcription factors

savePath = [writePath 'n008_' folder_prefix '02' filesep];
mkdir(savePath);
addpath(savePath)

[networkInfo, RSym8] = generate8StateNSNetworkBasic(networkInfoInit);

save([savePath 'networkInfo'],'networkInfo');    
networkInfo8 = networkInfo;
clear networkInfo
matlabFunction(RSym8,'File',[savePath 'RSymFun'],'Optimize',true,...
          'Vars',networkInfo8.sweepVarList);

%% iterate through higher-order models (3,4,&5 NS reactions)
% as we add additional reaction, this will multiply the state space by
% a factor of 2, since each state from smaller "base" network can coincide
% with the new factor being (1) engaged or (2) disengaged

n_ns_vec = 3:5;
n_state_vec = 2.^(1:4)*8;
networkInfoTemp = networkInfo8;

for n = 1%:length(n_state_vec)
    nStates = n_state_vec(n);
    % generate save path
    savePath = [writePath 'n' sprintf('%03d',nStates) '_' folder_prefix sprintf('%02d',n+2) filesep];
    mkdir(savePath);
    addpath(savePath)

    % call script to extend model
    [networkInfo, RSymTemp] = addNSReactionBasic(networkInfoTemp);
    save([savePath 'networkInfo'],'networkInfo');     
    networkInfoTemp = networkInfo;

    clear networkInfo
    matlabFunction(RSymTemp,'File',[savePath 'RSymFun'],'Optimize',true,...
              'Vars',networkInfoTemp.sweepVarList);
end        
