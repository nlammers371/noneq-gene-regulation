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
writePath = [DropboxFolder filesep 'sweeps03B_ir_vs_rate_cw_v2' filesep];
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
% NL: note that we're runnign for 100 iterations this time, since I saw
% that some eq sweeps were not converging
sweep_options = {'n_sim',50,'n_seeds',15,'n_iters_max',100, 'numerical_precision',10, ...
                'useParpool',1,'TauCycleTime',1,'downsample_output',ds_flag}; 
%%   
rate_index = find(strcmp(metric_names,'ProductionRate'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
precision_index = find(strcmp(metric_names,'Precision'));
phi_index = find(strcmp(metric_names,'Phi'));    
tau_index = find(strcmp(metric_names,'TauCycle'));  
cw_index = find(strcmp(metric_names,'CW'));  
ir_index = find(strcmp(metric_names,'IR'));        

% generate option vectors
% bs_vec = [1 2 3 4 5 1 1 1 1];
% lc_vec = [1 1 1 1 1 1 2 3 4];
% eq_vec = [1 1 1 1 1 0 0 0 0];
bs_vec = [1 2 3 4 1 1 1 1];
lc_vec = [1 1 1 1 1 2 3 4];
eq_vec = [1 1 1 1 0 0 0 0];
cw_vec = [10 1e2 1e3];
bs_vec_long = repelem(bs_vec,length(cw_vec));
lc_vec_long = repelem(lc_vec,length(cw_vec));
eq_vec_long = repelem(eq_vec,length(cw_vec));
cw_vec_long = repmat(cw_vec,1,length(bs_vec));


for i = 1:length(bs_vec_long)
    n_b = bs_vec_long(i);
    n_lc = lc_vec_long(i);
    eq_flag = eq_vec_long(i);
    cw = cw_vec_long(i);

    tic
    [sim_info, ~, sim_results] = ...
                        param_sweep_multi_v3([rate_index ir_index],...
                        functionPathCell{n_lc,n_b}, sweep_options{:},...
                        'equilibrium_flag', eq_flag, 'cw',cw);

%     if ds_flag
%         sim_results = sim_results_short;
%     end            
    suffix = '_neq';   
    if eq_flag
        suffix = '_eq';   
    end
    suffix = [suffix '_w' sprintf('%04d',cw)];
    % save
    disp('saving...')
    saveName = saveNameCell{n_lc,n_b};
    save([writePath 'sweep_info_' saveName suffix '.mat'],'sim_info')
    save([writePath 'sweep_results_' saveName suffix '.mat'],'sim_results', '-v7.3')

    toc;    

end    

