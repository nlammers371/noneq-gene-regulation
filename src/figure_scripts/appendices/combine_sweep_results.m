clear
close all

DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\';
DataPath1 = [DropboxFolder  'SweepOutput\sweeps03_info_vs_cw' filesep ];
DataPath2 = [DropboxFolder  'SweepOutput\sweeps03_info_vs_cw_rerun' filesep ];

DataPathOut = [DropboxFolder  'SweepOutput\sweeps03_info_vs_cw_cb' filesep ];
mkdir(DataPathOut)

% get lists of files
multi_bs_sweep_files1 = dir([DataPath1 'sweep_results*g01_cw1_eq*']);
multi_bs_info_files1 = dir([DataPath1 'sweep_info*g01_cw1_eq*']);

multi_lc_sweep_files1 = dir([DataPath1 'sweep_results*s01*_cw1.mat*']);
multi_lc_info_files1 = dir([DataPath1 'sweep_info*s01*_cw1.mat*']);

multi_bs_sweep_files2 = dir([DataPath2 'sweep_results*g01_cw1_eq*']);
multi_bs_info_files2 = dir([DataPath2 'sweep_info*g01_cw1_eq*']);

multi_lc_sweep_files2 = dir([DataPath2 'sweep_results*s01*_cw1_neq*']);
multi_lc_info_files2 = dir([DataPath2 'sweep_info*s01*_cw1_neq*']);

%% iterate through and concatenate
for i = 1:length(multi_bs_info_files1)
    % load BS files (1)
    load([DataPath1  multi_bs_info_files1(i).name])
    load([DataPath1  multi_bs_sweep_files1(i).name])
    
    sim_results_bs1 = sim_results;
    sim_info_bs1 = sim_info;
    
    clear sim_results
    clear sim_info
    % load BS files (2)
    load([DataPath2  multi_bs_info_files2(i).name])
    load([DataPath2  multi_bs_sweep_files2(i).name])
    
    sim_results_bs2 = sim_results;
    sim_info_bs2 = sim_info;
    
    % combine
    sim_info = sim_info_bs1;
    sim_info.n_iter_vec = [sim_info_bs1.n_iter_vec sim_info_bs2.n_iter_vec];
    sim_info.convergence_flag_vec = [sim_info_bs1.convergence_flag_vec sim_info_bs2.convergence_flag_vec];
    sim_info.n_sim = length(sim_info.n_iter_vec);
    
    sim_results_bs2.areaIDVec = sim_results_bs2.areaIDVec + sim_info_bs1.n_sim;
    sim_results = struct;
    sim_results.metric_array = vertcat(sim_results_bs1.metric_array, sim_results_bs2.metric_array);
    sim_results.rate_array = vertcat(sim_results_bs1.rate_array, sim_results_bs2.rate_array);
    
    % save
    save([DataPathOut multi_bs_info_files2(i).name],'sim_info')
    save([DataPathOut multi_bs_sweep_files2(i).name],'sim_results')
    
end  

%%

for i = 1:4%length(multi_lc_info_files1)
    % load BS files (1)
    load([DataPath1  multi_lc_info_files1(i).name])
    load([DataPath1  multi_lc_sweep_files1(i).name])
    
    sim_results_lc1 = sim_results;
    sim_info_lc1 = sim_info;
    
    clear sim_results
    clear sim_info
    % load BS files (2)
    load([DataPath2  multi_lc_info_files2(i).name])
    load([DataPath2  multi_lc_sweep_files2(i).name])
    
    sim_results_lc2 = sim_results;
    sim_info_lc2 = sim_info;
    
    % combine
    sim_info = sim_info_lc1;
    sim_info.n_iter_vec = [sim_info_lc1.n_iter_vec sim_info_lc2.n_iter_vec];
    sim_info.convergence_flag_vec = [sim_info_lc1.convergence_flag_vec sim_info_lc2.convergence_flag_vec];
    sim_info.n_sim = length(sim_info.n_iter_vec);
    
    sim_results_lc2.areaIDVec = sim_results_lc2.areaIDVec + sim_info_lc1.n_sim;
    sim_results = struct;
    sim_results.metric_array = vertcat(sim_results_lc1.metric_array, sim_results_lc2.metric_array);
    sim_results.rate_array = vertcat(sim_results_lc1.rate_array, sim_results_lc2.rate_array);
    
    % save
    save([DataPathOut multi_lc_info_files2(i).name],'sim_info')
    save([DataPathOut multi_lc_sweep_files2(i).name],'sim_results')
    
end  