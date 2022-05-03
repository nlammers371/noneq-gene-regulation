% script to call core parameter sweep function to examine tradeoffs between
% different transcritional input/output behaviors

clear 
close all
addpath(genpath('../utilities/'))

% as we add additional reaction, this will multiply the state space by
% a factor of 2, since each state from smaller "base" network can coincide
% with the new factor being (1) engaged or (2) disengaged

% set basic parameters
n_bs_vec = 1:5;
n_g_vec = 1:5;
ns_flag = 0;

% Set Dropbox directory
DropboxFolder = 'S:\Nick\Dropbox\Nonequilibrium\Nick\SweepOutput';
writePath = [DropboxFolder filesep 'sweeps01_info_vs_energy_v3' filesep];
mkdir(writePath);

% this contains paths used to address correct functions
functionPathCell = cell(length(n_g_vec),length(n_bs_vec));
saveNameCell = cell(length(n_g_vec),length(n_bs_vec));

% get metric names
[~,~,metric_names] = calculateMetricsNumeric_v3([]);

% generate path to save metric functions 
subfolderName = 'numeric';
sourcePath = handlePathOptions(subfolderName);


for m = 1:length(n_g_vec)
    for n = 1:length(n_bs_vec) 

        % add path to functions        
        readName = ['s' sprintf('%02d',n_bs_vec(n)) '_ns00_g'...
                  sprintf('%02d',n_g_vec(m)) '_cw' num2str(ns_flag)];
        functionPath = [sourcePath readName filesep];
        functionPathCell{m,n} = functionPath;
        saveNameCell{m,n} = readName;

    end
end    
                      
% set sim options
sweep_options = {'n_sim',200,'n_seeds',15,'n_iters_max',50, 'numerical_precision',10, 'useParpool',1};
%%   
phi_index = find(strcmp(metric_names,'Phi'));    
ir_index = find(strcmp(metric_names,'IR'));
        

for m = 1:length(n_g_vec)
    if m < 2
        bs_list = 1:5;
    else
        bs_list = 1;%:2;
    end
    for n = bs_list
                              
        % call sweep function
        tic
        [sim_info, ~, sim_results] = ...
                            param_sweep_multi_v3([phi_index ir_index],...
                            functionPathCell{m,n}, sweep_options{:},...
                            'equilibrium_flag',false,'TauCycleTime',1,...
                            'paramBounds',paramBounds,'downsample_output',true);
                          
        % save
        disp('saving...')
        saveName = saveNameCell{m,n};
        save([writePath 'sweep_info_' saveName '.mat'],'sim_info')
        save([writePath 'sweep_results_' saveName '.mat'],'sim_results', '-v7.3')
        
        toc;    

    end
end    
