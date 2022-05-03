% This generates models with rate constraints intended to reflect realistic
% constraints for different preposed regulatory mechanisms
clear
close all

% generate baseline 6 state network, which corresponds to system with 1
% binding site that is subjected to both correct and incorrect activator
% binding
[~, activeStatesBase, n_wrong_bound_init, n_right_bound_init, activity_vec_init,...
          transitionInfoInit,permittedConnectionsInit]...
            = generateBasicSixStateNetwork;          

%% generate "Nucleosome" model
% This corresponds to network with two binding sites, one of which is
% assumed to be non-specific
savePath = '../utilities/metricFunctions/n6_OR_Nucleosome_NUM/';
mkdir(savePath);
addpath(savePath)

[networkInfo, RSymNuc] = generate6StateNucleosomeNetwork(n_wrong_bound_init,...
                                n_right_bound_init,activity_vec_init,...
                                transitionInfoInit, savePath);

save([savePath 'networkInfo'],'networkInfo');    
matlabFunction(RSymNuc,'File',[savePath 'RSymFun'],'Optimize',true,...
          'Vars',networkInfo.sweepVarList);
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

    [networkInfo, RSymTemp] = addSpecBindingSite_v2(networkInfoTemp);
    save([savePathTemp 'networkInfo'],'networkInfo');     
    networkInfoTemp = networkInfo;

    clear networkInfo
    matlabFunction(RSymTemp,'File',[savePathTemp 'RSymFun'],'Optimize',true,...
              'Vars',networkInfoTemp.sweepVarList);
end        
