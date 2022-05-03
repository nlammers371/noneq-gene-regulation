% script to examine sharpness-precision tradeoffs in systems without
% activator-activator cooperativity



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
writePath = [DropboxFolder filesep 'sweeps02_sharpness_vs_precision_no_pairwise' filesep];
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
sweep_options = {'n_sim',24,'n_seeds',15,'n_iters_max',50, 'numerical_precision',10, ...
            'useParpool',1,'downsample_output',ds_flag==1,'TauCycleTime',1};
%%   
rate_index = find(strcmp(metric_names,'ProductionRate'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
precision_index = find(strcmp(metric_names,'Precision'));
phi_index = find(strcmp(metric_names,'Phi'));    
tau_index = find(strcmp(metric_names,'TauCycle'));    
ir_index = find(strcmp(metric_names,'IR'));

[sim_info, ~, ~] = param_sweep_multi_v3([sharpness_index precision_index],...
                                functionPathCell{m,n}, 'n_sim',0);
                              
%%                              
       
% paramBounds = repmat([-5; 5],1,10);
% paramBounds(1,6) = 0;
% paramBounds(2,7) = 0;

% results_struct = struct;

for m = 1%length(n_g_vec)

    for n = 1:5
                       
        [sim_info, ~, ~] = param_sweep_multi_v3([sharpness_index precision_index],...
                                functionPathCell{m,n}, 'n_sim',0);
      
        sweepVarStrings = sim_info.sweepVarStrings;
        sweepFlags = true(size(sweepVarStrings));
        cr_index = find(strcmp(sweepVarStrings,'cr'));
        wmp_index = find(strcmp(sweepVarStrings,'wmp'));
        sweepFlags([wmp_index cr_index ]) = false;                
                              
        % call sweep function
        tic
        if ~eq_only_flag
            [sim_info, sim_results_long, sim_results] = ...
                                param_sweep_multi_v3([sharpness_index precision_index],...
                                functionPathCell{m,n}, sweep_options{:},'equilibrium_flag',false,'sweepFlags',sweepFlags);
            if ~ds_flag
                sim_results = sim_results_long;
            end
%             save
            disp('saving neq results...')
            saveName = saveNameCell{m,n};
            save([writePath 'sweep_info_' saveName '_neq.mat'],'sim_info');
            save([writePath 'sweep_results_' saveName '_neq.mat'],'sim_results', '-v7.3');
        end
        if ~neq_only_flag
            [sim_info, sim_results_long, sim_results] = ...
                                param_sweep_multi_v3([sharpness_index precision_index],...
                                functionPathCell{m,n}, sweep_options{:},'paramBounds',paramBounds,'equilibrium_flag',true,'sweepFlags',sweepFlags);
            if ~ds_flag
                sim_results = sim_results_long;
            end
            % save
            disp('saving eq results...')
            saveName = saveNameCell{m,n};
            save([writePath 'sweep_info_' saveName '_eq.mat'],'sim_info');
            save([writePath 'sweep_results_' saveName '_eq.mat'],'sim_results', '-v7.3');
        end          
        toc;    
    end
end    

