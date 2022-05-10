% script to call core parameter sweep function to examine tradeoffs between
% different transcritional input/output behaviors

clear 
close all
addpath(genpath('../utilities/'))

% set basic parameters
n_bs_vec = 1:5;
n_g_vec = 1:5;
ns_flag = 1;
ds_flag = 1;

% Set Dropbox directory
DropboxFolder = 'S:\Nick\Dropbox\Nonequilibrium\Nick\SweepOutput';
writePath = [DropboxFolder filesep 'sweeps03_info_vs_cw_rerun' filesep];
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
sweep_options = {'n_sim',1e3,'n_seeds',15,'n_iters_max',100, 'numerical_precision',10, ...
                'useParpool',1,'TauCycleTime',1,'downsample_output',ds_flag};
%%   
rate_index = find(strcmp(metric_names,'ProductionRate'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
precision_index = find(strcmp(metric_names,'Precision'));
phi_index = find(strcmp(metric_names,'Phi'));    
tau_index = find(strcmp(metric_names,'TauCycle'));  
cw_index = find(strcmp(metric_names,'CW'));  
ir_index = find(strcmp(metric_names,'IR'));
        

for equilibrium_flag = 0
    for m = 1:4
%         if false%m == 1
%             bs_list = 1:5;
%         else
%             bs_list = 2:5;
%         end
        for n = 1
            % call sweep function
            tic
            [sim_info, sim_results, sim_results_short] = ...
                                param_sweep_multi_v3([cw_index ir_index],...
                                functionPathCell{m,n}, sweep_options{:},...
                                'equilibrium_flag',equilibrium_flag);

            if ds_flag
                sim_results = sim_results_short;
            end            
            suffix = '_neq';
            if equilibrium_flag
                suffix = '_eq';
            end
            
            % save
            disp('saving...')
            saveName = saveNameCell{m,n};
            save([writePath 'sweep_info_' saveName suffix '.mat'],'sim_info')
            save([writePath 'sweep_results_' saveName suffix '.mat'],'sim_results', '-v7.3')

            toc;    
    
        end
    end    
end
