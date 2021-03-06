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
ds_flag = 1;
eq_only_flag = 0;
neq_only_flag = 0;
% Set Dropbox directory
DropboxFolder = 'S:\Nick\Dropbox\Nonequilibrium\Nick\SweepOutput';
writePath = [DropboxFolder filesep 'sweeps02D_S_vs_r_prec12' filesep];
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
sweep_options = {'n_sim',50,'n_seeds',15,'n_iters_max',200, 'numerical_precision',12, ...
            'useParpool',1,'downsample_output',ds_flag,'TauCycleTime',1};
%%   
rate_index = find(strcmp(metric_names,'ProductionRate'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));

% cw_val = 1000;

for m = 1%[4:-1:1]%1:length(n_g_vec)-1
%     if m == 1
%         bs_list = 1; 
%     else
%         bs_list = 1;
%     end
    for n = 5:-1:1%bs_list                              
        % call sweep function
        tic

        [sim_info, sim_results_long, sim_results] = ...
                            param_sweep_multi_v3([rate_index sharpness_index],...
                            functionPathCell{m,n}, sweep_options{:},'equilibrium_flag', true);
        if ~ds_flag
            sim_results = sim_results_long;
        end
        % save
        disp('saving neq results...')
        saveName = saveNameCell{m,n};
        save([writePath 'sweep_info_' saveName '_eq.mat'],'sim_info');
        save([writePath 'sweep_results_' saveName '_eq.mat'],'sim_results', '-v7.3');
      
        toc        
    end
end    
