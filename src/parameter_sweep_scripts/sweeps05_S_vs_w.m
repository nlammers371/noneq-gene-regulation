% Sweep observed sharpness as a function of w/c

clear 
close all
addpath(genpath('../utilities/'))


% set basic parameters
n_bs_vec = 1:5;
n_g_vec = 1:5;
ns_flag = 1;
ds_flag = 1;
% hm_flag = 1;
% Set Dropbox directory
DropboxFolder = 'S:\Nick\Dropbox\Nonequilibrium\Nick\SweepOutput';
writePath = [DropboxFolder filesep 'sweeps05C_S_vs_w_rerun' filesep];
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
sweep_options = {'n_sim',500,'n_seeds',15,'n_iters_max',50, 'numerical_precision',10, 'useParpool',1,'ds_factor',100};

%%   
cw_index = find(strcmp(metric_names,'CW'));  
ir_index = find(strcmp(metric_names,'IR'));

for equilibrium_flag = 0
    for m = fliplr(1:4)
        if m == 1
            bs_list = fliplr(1:5);
        else
            bs_list = 1;
        end
        
        for n = 1
            % call sweep function
            tic
            [sim_info, sim_results, sim_results_short] = ...
                                param_sweep_multi_v3([cw_index ir_index],...
                                functionPathCell{m,n}, sweep_options{:},...
                                'equilibrium_flag',equilibrium_flag,'TauCycleTime',1,...
                                'downsample_output',ds_flag);

            if ds_flag
                sim_results = sim_results_short;
            end            
            suffix = '';
            if equilibrium_flag
                suffix = '_eq';
            end
            % save
            disp('saving...')
            saveName = saveNameCell{m,n};
            save([writePath 'sweep_info_' saveName suffix '.mat'],'sim_info')
            save([writePath 'sweep_results_' saveName suffix '.mat'],'sim_results', '-v7.3')

            toc;    
    %         iter = iter + 1;
        end
    end    
end
