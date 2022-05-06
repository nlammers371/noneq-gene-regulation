% Plot results for IR vs energy for higher order models
clear 
close all
addpath(genpath('../utilities/'))

% %%%%%%%%%%%%%%%%  set relative read and write paths %%%%%%%%%%%%%%%%%%%%
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\';
DataPath = [DropboxFolder  'SweepOutput\sweeps03B_ir_vs_rate_cw' filesep ];
DataPath_f5 = [DropboxFolder  'SweepOutput\sweeps03_info_vs_cw' filesep ];
FigPath = [DropboxFolder '\manuscript\appendices' filesep 'sweep_algorithm' filesep];
mkdir(FigPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% load numerical sweep results sets MULTI BS
[~,~,metric_names_num] = calculateMetricsNumeric_v3([]);
n_ir = 100;
ir_index = find(strcmp(metric_names_num,'IR'));
w_index = find(strcmp(metric_names_num,'CW'));
%%%%%%%%%
% Multi BS
%%%%%%%%%

% get list of sweep results files with only 1 genera TF reaction
ref_sweep_files = dir([DataPath 'sweep_results*']);
ref_info_files = dir([DataPath 'sweep_info*']);

ref_struct = struct; 
% load
for f = 1:length(ref_info_files)
  
    % load files
    load([DataPath  ref_info_files(f).name])
    load([DataPath  ref_sweep_files(f).name])
    
    % add main structure
    ref_struct(f).sweep_results = sim_results;
    ref_struct(f).sweep_info = sim_info;        
    
    % add key info
    ref_struct(f).eq_flag = sim_info.equilibrium_flag;
    ref_struct(f).w = sim_info.cw;
    
    s_ind = strfind(ref_sweep_files(f).name,'_s');
    ref_struct(f).nb = str2double(ref_sweep_files(f).name(s_ind+2:s_ind+3));
    
    g_ind = strfind(ref_sweep_files(f).name,'_g');
    ref_struct(f).nlc = str2double(ref_sweep_files(f).name(g_ind+2:g_ind+3));
    
    % calculate max info
    ref_struct(f).max_ir = nanmax(sim_results.metric_array(:,ir_index))*log2(exp(1));
    
    clear sim_info
    clear sim_results
end            

% generate helper vectors
w_val_vec = [ref_struct.w];
nb_val_vec = [ref_struct.nb];
nlc_val_vec = [ref_struct.nlc];
eq_flag_vec = [ref_struct.eq_flag];
max_ir_vec = [ref_struct.max_ir];

%% Load IR vs. W results used in main text Figure 5
w_index_vec = unique(w_val_vec);

multi_bs_sweep_files_f5 = dir([DataPath_f5 'sweep_results*g01*eq*']);
multi_bs_info_files_f5 = dir([DataPath_f5 'sweep_info*g01*eq*']);

multi_lc_sweep_files_f5 = dir([DataPath_f5 'sweep_results*s01*cw1.mat']);
multi_lc_info_files_f5 = dir([DataPath_f5 'sweep_info*s01*cw1.mat']);

master_struct_multi_lc = struct;
master_struct_multi_bs = struct;

for f = 1:length(multi_bs_sweep_files_f5)
  
    % load files
    load([DataPath_f5  multi_bs_info_files_f5(f).name])
    load([DataPath_f5  multi_bs_sweep_files_f5(f).name])
    
    master_struct_multi_bs(f).sweep_results_f5 = sim_results;
    master_struct_multi_bs(f).sweep_info_f5 = sim_info;        
    
    clear sim_results
    clear sim_info
    
    % calculate max IR values for each relevant W level
    w_vec = 10.^master_struct_multi_bs(f).sweep_results_f5.metric_array(:,w_index);
    ir_vec = master_struct_multi_bs(f).sweep_results_f5.metric_array(:,ir_index)*log2(exp(1));
    
    master_struct_multi_bs(f).nb_vec = repelem(f,length(w_index_vec));
    master_struct_multi_bs(f).nlc_vec = repelem(1,length(w_index_vec));
    master_struct_multi_bs(f).eq_flag_vec = repelem(master_struct_multi_bs(f).sweep_info_f5.equilibrium_flag,length(w_index_vec));
    ir_max_vec = NaN(size(w_index_vec));
    master_struct_multi_bs(f).w_ref = w_index_vec;
    for i = 1:length(w_index_vec)
        w_filter = w_vec>=w_index_vec(i);
        ir_max_vec(i) = nanmax(ir_vec(w_filter));
    end
    master_struct_multi_bs(f).ir_max_vec = ir_max_vec;
end            

%%%%%%%%%
% Multi-LC
%%%%%%%%%

for f = 1:4
  
    load([DataPath_f5 multi_lc_info_files_f5(f).name])
    load([DataPath_f5  multi_lc_sweep_files_f5(f).name])
    
    master_struct_multi_lc(f).sweep_results_f5 = sim_results;
    master_struct_multi_lc(f).sweep_info_f5 = sim_info;
    
    clear sim_results
    clear sim_info
    
    % calculate max IR values for each relevant W level
    w_vec = 10.^master_struct_multi_lc(f).sweep_results_f5.metric_array(:,w_index);
    ir_vec = master_struct_multi_lc(f).sweep_results_f5.metric_array(:,ir_index)*log2(exp(1));
    
    master_struct_multi_lc(f).nb_vec = repelem(1,length(w_index_vec));
    master_struct_multi_lc(f).nlc_vec = repelem(f,length(w_index_vec));
    master_struct_multi_lc(f).eq_flag_vec = repelem(master_struct_multi_lc(f).sweep_info_f5.equilibrium_flag,length(w_index_vec));
    ir_max_vec = NaN(size(w_index_vec));
    master_struct_multi_lc(f).w_ref = w_index_vec;
    for i = 1:length(w_index_vec)
        w_filter = w_vec>=w_index_vec(i);
        ir_max_vec(i) = nanmax(ir_vec(w_filter));
    end
    master_struct_multi_lc(f).ir_max_vec = ir_max_vec;
end 

%% Identify matches and compare IR values

% generate helper vecs
max_ir_vec_full = [[master_struct_multi_lc.ir_max_vec] [master_struct_multi_bs.ir_max_vec]];
w_vec_full = [[master_struct_multi_lc.w_ref] [master_struct_multi_bs.w_ref]];
eq_flag_vec_full = [[master_struct_multi_lc.eq_flag_vec] [master_struct_multi_bs.eq_flag_vec]];
nb_vec_full = [[master_struct_multi_lc.nb_vec] [master_struct_multi_bs.nb_vec]];
nlc_vec_full = [[master_struct_multi_lc.nlc_vec] [master_struct_multi_bs.nlc_vec]];


full_sweep_ir_vals = NaN(size(max_ir_vec));
for i = 1:length(max_ir_vec)
    match_filter = w_vec_full==w_val_vec(i) & eq_flag_vec_full==eq_flag_vec(i) &...
                   nb_vec_full==nb_val_vec(i) & nlc_vec_full==nlc_val_vec(i);
                 
    if sum(match_filter) > 1
        error('too many matches')
    elseif sum(match_filter) == 1
        full_sweep_ir_vals(i) = max_ir_vec_full(match_filter);
    end
end    
