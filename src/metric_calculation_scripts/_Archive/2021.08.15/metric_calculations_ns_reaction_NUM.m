% This generates increasingly complex architecture files for use in
% numerical parameter sweeps. We assume 1 specific binding site and an
% increasing number of distinct, non-specific reactions 
clear
close all

% generate baseline 6 state network, which corresponds to system with 1
% binding site that is subjected to both correct and incorrect activator
% binding
[~, activeStatesBase, n_wrong_bound_init, n_right_bound_init, activity_vec_init,...
          transitionInfoInit,permittedConnectionsInit]...
            = generateBasicSixStateNetwork;
          
% excise non-specific states
baseNum = 4;          
n_right_bound_init = n_right_bound_init(1:baseNum);
activity_vec_init = activity_vec_init(1:baseNum);
transitionInfoInit = transitionInfoInit(1:baseNum,1:baseNum);
permittedConnectionsInit = permittedConnectionsInit(1:baseNum,1:baseNum);

%% generate full 8 network info
% This corresponds to network with two binding sites, one of which is
% assumed to be non-specific
savePath = '../utilities/metricFunctions/n8_OR_NS_NUM/';
mkdir(savePath);
addpath(savePath)

[networkInfo, RSym8] = generate8StateNSNetwork(baseNum...
                                ,n_right_bound_init,activity_vec_init,...
                                transitionInfoInit,permittedConnectionsInit);

save([savePath 'networkInfo'],'networkInfo');    
networkInfo8 = networkInfo;
clear networkInfo
matlabFunction(RSym8,'File',[savePath 'RSymFun'],'Optimize',true,...
          'Vars',networkInfo8.sweepVarList);
% copyfile(fourStateRefPath, [savePath 'n4_OR'])
%% iterate through higher order models (3,4,5,&6 sites)
% as we add additional binding sites, this will multiply the state space by
% a factor of 2, since each state from smaller "base" network can coincide
% with the new site being (1) empty or (2) bound by the correct factor

n_site_vec = 3:6;
n_state_vec = 2.^(1:4)*8;
networkInfoTemp = networkInfo8;
for n = 1:length(n_state_vec)
    nStates = n_state_vec(n);
    % generate save path
    savePathTemp = ['../utilities/metricFunctions/n' num2str(nStates) '_OR_SPEC_NUM/'];
    mkdir(savePathTemp);
    addpath(savePathTemp)

    [networkInfo, RSymTemp] = addSpecBindingSite(networkInfoTemp);
    save([savePathTemp 'networkInfo'],'networkInfo');     
    networkInfoTemp = networkInfo;

    clear networkInfo
    matlabFunction(RSymTemp,'File',[savePathTemp 'RSymFun'],'Optimize',true,...
              'Vars',networkInfoTemp.sweepVarList);
end        
