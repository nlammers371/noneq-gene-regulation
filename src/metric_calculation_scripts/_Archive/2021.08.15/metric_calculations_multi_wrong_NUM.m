% This script explores possibility of keeping SS vec in symbolic form
clear
close all

% generate baseline 6 state network, which corresponds to system with 1
% binding site that is subjected to both correct and incorrect activator
% binding
[baseNum, activeStatesBase, n_wrong_bound_init,n_right_bound_init,activity_vec_init,...
          transitionInfoInit,permittedConnectionsInit]...
            = generateBasicSixStateNetwork;

% generate full 18 network info
% This corresponds to network with two binding sites, one of which is
% assumed to be non-specific
savePath = '../utilities/metricFunctions/n18_OR_NUM/';
mkdir(savePath);
addpath(savePath)

[networkInfo, RSym18] = generate18StateNetwork_v3(baseNum, n_wrong_bound_init...
                                  ,n_right_bound_init,activity_vec_init,...
                                   transitionInfoInit,permittedConnectionsInit);

save([savePath 'networkInfo'],'networkInfo');    
networkInfo18 = networkInfo;
clear networkInfo
matlabFunction(RSym18,'File',[savePath 'RSymFun'],'Optimize',true,...
          'Vars',networkInfo18.sweepVarList);
        
% generate 6 state network using 18 state base
savePath = '../utilities/metricFunctions/n6_OR_NUM/';
mkdir(savePath);
addpath(savePath)

networkInfo = truncate18To6(networkInfo18);

save([savePath 'networkInfo'],'networkInfo');    
networkInfo6 = networkInfo;
clear networkInfo
matlabFunction(networkInfo6.RSym,'File',[savePath 'RSymFun'],'Optimize',true,...
          'Vars',networkInfo6.sweepVarList);

% generate 3x18=54 state network 
% as we add additional binding sites, this will multiply the state space by
% a factor of 3, since each state from smaller "base" network can coincide
% with the new site being (1) empty, (2) bound by the correct factor, (3) bound
% by the incorrect factor
savePath54 = '../utilities/metricFunctions/n54_OR_NUM/';
mkdir(savePath54);
addpath(savePath54)

[networkInfo, RSym54] = addBindingSite_v3(networkInfo18);
save([savePath54 'networkInfo'],'networkInfo');     
networkInfo54 = networkInfo;

clear networkInfo
matlabFunction(RSym54,'File',[savePath54 'RSymFun'],'Optimize',true,...
          'Vars',networkInfo54.sweepVarList);
        
% generate 3x54=162 state network (4 binding sites)

savePath162 = '../utilities/metricFunctions/n162_OR_NUM/';
mkdir(savePath162);
addpath(savePath162)

[networkInfo, RSym162] = addBindingSite_v3(networkInfo54);
save([savePath162 'networkInfo'],'networkInfo');     
networkInfo162 = networkInfo;
clear networkInfo
matlabFunction(RSym162,'File',[savePath162 'RSymFun'],'Optimize',true,...
          'Vars',networkInfo162.sweepVarList);        