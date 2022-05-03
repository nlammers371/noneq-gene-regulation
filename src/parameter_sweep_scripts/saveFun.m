function saveFun(suffix, writePath, sweep_info_neq, sweep_results_neq, sweep_info_eq, sweep_results_eq, sweep_info_eq0, sweep_results_eq0)

    save([writePath 'sweep_info_neq_' suffix '.mat'], 'sweep_info_neq')
    save([writePath 'sweep_results_neq_' suffix '.mat'], 'sweep_results_neq')
    
    save([writePath 'sweep_info_eq_' suffix '.mat'], 'sweep_info_eq')
    save([writePath 'sweep_results_eq_' suffix '.mat'], 'sweep_results_eq')
    
    save([writePath 'sweep_info_eq0_' suffix '.mat'], 'sweep_info_eq0')
    save([writePath 'sweep_results_eq0_' suffix '.mat'], 'sweep_results_eq0')