% This generates increasingly complex architecture files for use in
% numerical parameter sweeps. In this case, we allow general factors to be
% non-identical

clear
close all

% generate path to save metric functions 
subfolderNameFull = 'numeric/unique/';
writePathFull = handlePathOptions(subfolderNameFull);

subfolderNameG = 'numeric/unique_g/';
writePathG = handlePathOptions(subfolderNameG);

subfolderNameA = 'numeric/unique_a/';
writePathA = handlePathOptions(subfolderNameA);

% first, generate baseline 6 state network that serves as kernel for all
% architectures
networkInfoInit = generateBaselineSixStateNetwork;

% generate generic prefix for function directories 
% folder_string = 'ns00_g01';       
%% generate systems with 2-10 binding sites
% as we add additional reaction, this will multiply the state space by
% a factor of 2, since each state from smaller "base" network can coincide
% with the new factor being (1) engaged or (2) disengaged

n_bs_vec = 1:2;
n_g_vec = 1:2;
% n_bs_cw = [6 12 20 30 42]/2;
% n_state_vec = (n_bs_vec+1).*(n_g_vec+1)';
% n_state_vec_cw = (n_bs_cw).*(n_g_vec+1)';

for ns_flag = 1%0:1
    for m = 1:length(n_g_vec)
        for n = 1:length(n_bs_vec)                                  

%             % call script to genereate model
%             
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             % Start with fully general
%             if ns_flag == 0
%                 error('not yet supported');%[networkInfo, RSym, RSymFull] = makeActivatorGeneralNetwork(networkInfoInit,n_bs_vec(n),n_g_vec(m));
%             else
%                 [networkInfo, RSym, RSymFull] = makeActivatorGeneralNetworkNSDiffAll(networkInfoInit,n_bs_vec(n),n_g_vec(m));
%             end
% 
%             nStates = size(RSym,1);
%             
%             savePathFull = [writePathFull 's' sprintf('%02d',n_bs_vec(n)) '_ns00' ...
%                         '_g' sprintf('%02d',n_g_vec(m)) '_cw' num2str(ns_flag) filesep];
%                       
%             mkdir(savePathFull);
%             addpath(savePathFull)
%             
%             save([savePathFull 'networkInfo'],'networkInfo');    
%                     
%             matlabFunction(RSym,'File',[savePathFull 'RSymFun'],'Optimize',true,...
%                                         'Vars',networkInfo.sweepVarList);             
% 
%             
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             % Generic activator binding but unique G reactions
%             if ns_flag == 0
% %                 [networkInfo, RSym, RSymFull] = makeActivatorGeneralNetwork(networkInfoInit,n_bs_vec(n),n_g_vec(m));
%             else
%                 [networkInfo, RSym, RSymFull] = makeActivatorGeneralNetworkNSDiffGen(networkInfoInit,n_bs_vec(n),n_g_vec(m));
%             end
% 
%             nStates = size(RSym,1);
%             
%             savePathG = [writePathG 's' sprintf('%02d',n_bs_vec(n)) '_ns00' ...
%                         '_g' sprintf('%02d',n_g_vec(m)) '_cw' num2str(ns_flag) filesep];
%                       
%             mkdir(savePathG);
%             addpath(savePathG)
%             
%             save([savePathG 'networkInfo'],'networkInfo');    
%                     
%             matlabFunction(RSym,'File',[savePathG 'RSymFun'],'Optimize',true,...
%                                         'Vars',networkInfo.sweepVarList); 
%                                       
%                                       
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Generic general tf binding but unique A reactions
            if ns_flag == 0
%                 [networkInfo, RSym, RSymFull] = makeActivatorGeneralNetwork(networkInfoInit,n_bs_vec(n),n_g_vec(m));
            else
                [networkInfo, RSym, RSymFull] = makeActivatorGeneralNetworkNSDiffAct(networkInfoInit,n_bs_vec(n),n_g_vec(m));
            end

            nStates = size(RSym,1);
            
            savePathA = [writePathA 's' sprintf('%02d',n_bs_vec(n)) '_ns00' ...
                        '_g' sprintf('%02d',n_g_vec(m)) '_cw' num2str(ns_flag) filesep];
                      
            mkdir(savePathA);
            addpath(savePathA)
            
            save([savePathA 'networkInfo'],'networkInfo');    
                    
            matlabFunction(RSym,'File',[savePathA 'RSymFun'],'Optimize',true,...
                                        'Vars',networkInfo.sweepVarList);

        end        
    end
end    