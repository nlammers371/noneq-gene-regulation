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
% n_bs_cw = [6 12 20 30 42]/2;

for ns_flag = 1%0:1
      
    for m = 2%:5%m_max
%         if m < 3
%             bs_list = 1:5;
%         else
%             bs_list = 1:2;
%         end
        bs_list = 5;%3:4;
        for n = bs_list
            
            savePath = [writePath 's' sprintf('%02d',n_bs_vec(n)) '_ns00' ...
                        '_g' sprintf('%02d',n_g_vec(m)) '_cw' num2str(ns_flag) filesep];
            mkdir(savePath);
            addpath(savePath)

            % call script to extend model
            if ns_flag == 0
                [networkInfo, RSym, RSymFull] = makeActivatorGeneralNetwork(networkInfoInit,n_bs_vec(n),n_g_vec(m));
            else
                [networkInfo, RSym, RSymFull] = makeActivatorGeneralNetworkNS(networkInfoInit,n_bs_vec(n),n_g_vec(m));
            end

            save([savePath 'networkInfo'],'networkInfo');     
        %     clear networkInfo
            matlabFunction(RSym,'File',[savePath 'RSymFun'],'Optimize',true,...
                                        'Vars',networkInfo.sweepVarList);             

            if size(RSym,1) <= 16
                matlabFunction(RSymFull,'File',[savePath 'RSymFunFull'],'Optimize',false,...
                                      'Vars',networkInfo.sweepVarList);
            end

        end        
    end
end    