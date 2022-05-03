% script to estimate decision time limits for different organisms using
% outside data sources
clear 
close all
addpath(genpath('../utilities/'))

% basic stats
name_cell = {'Drosophila (early development)', 'C. elegans (adult)', 'Arabidopsis', 'Human (adult)', 'Mouse (adult)'};

%%%%%%%%%%%%%%%
% mRNA decay constants (in minutes)
%%%%%%%%%%%%%%%
% Note 1: upper bound for dmel is taken from length of nc14. Lower bound is
        % taken from (approximate) timescale of mRNA half-life (see Appendix B of
        % Lammers et all 2020 PNAS)
% Note 2: C elegans ub reflect cell cycle time estimates cited in: Hubbard, E. J. A. The C. elegans germ line: a model for stem cell biology. Dev. Dyn. 236, 3343 (2007).
% Note 3: C elegans lb reflect HL data for rpl-7A gene taken from panel k of Figure 4 in: Son, H. G. et al. RNA surveillance via nonsense-mediated mRNA decay is crucial for longevity in daf-2/insulin/IGF-1 mutant C. elegans. Nat. Commun. 2017 81 8, 1–11 (2017).
% Note 4: Bounds for Arabidopsis are sited as "minutes to days"; however
        % cell cylcle is given as 5.9 hours. We take this as upper bound and set
        % lower bound to 30 min
% Note 5: Remaining bounds are taken from HL estimates provided in Table 1 of: 1. Pérez-Ortín, J. E., Alepuz, P., Chávez, S. & Choder, M. Eukaryotic mRNA Decay: Methodologies, Pathways, and Links to Other Stages of Gene Expression. (2013). doi:10.1016/j.jmb.2013.02.029        

mRNA_half_life_ub = [60*log(2), 24*60*log(2), 5.9*60, 48*60, 30*60]/log(2);
mRNA_half_life_lb = [7*log(2),  6*60,         30,     4*60, 30]/log(2);

%%%%%%%%%%%%%%%
% Burst cycle time ranges 
%%%%%%%%%%%%%%%
% Note 1: Arabidopsis timescale was estimated by looking at raw data from: 1. Alamos, S., Reimer, A., Niyogi, K. K. & Garcia, H. G. Quantitative imaging of RNA polymerase II activity in plants reveals the single-cell basis of tissue-wide transcriptional dynamics. Nat. Plants 2021 78 7, 1037–1049 (2021).
% Note 2: C. elegans cycle times are taken from: Lee, C. H., Shin, H. & Kimble, J. Dynamics of Notch-Dependent Transcriptional Bursting in Its Native Context. Dev. Cell 50, 426-435.e4 (2019).
% Note 3: All other times are taken from Table 1 of Appendix A of: https://doi.org/10.1016/j.ceb.2020.08.001
burst_cycle_time_ub = [10, 105.3, 60, 3*60, 3*60];
burst_cycle_time_lb = [2,  60.5,  30, 60,   30];

%%%%%%%%%%%%%%%
% Number of TF species
%%%%%%%%%%%%%%%
% Note 1: dmel estimate comes from counting number of distinct TF species
        % included in SiteOut binding site removal software: https://doi.org/10.1371/journal.pone.0151740
% Note 2: Remainder come from database fo TF abundances accompanying this work: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2995046/  
% Note 3: Multiplicative factors are applied to account for variability in
        % TF expression levels across different cells. They are intended
        % for illustrative purposes and do not represent a quantitative
        % eistmate of TF expression variance
       
n_tfs_ub = [47, 698, 1356, 1508, 1426]*sqrt(10);
n_tfs_lb = [47, 698, 1356, 1508, 1426]/sqrt(10);

% generate estimates and (rough) uncertainty ranges
decision_limit_info = struct;

% cw
decision_limit_info.cw_ub = n_tfs_ub;
decision_limit_info.cw_lb = n_tfs_lb;

% n burst cycles for decision
decision_limit_info.n_cycles_ub = mRNA_half_life_ub./burst_cycle_time_lb;
decision_limit_info.n_cycles_lb = mRNA_half_life_lb./burst_cycle_time_ub;

save('decision_limit_info.mat','decision_limit_info')