function metric = rate_opt_general(opt_rates,opt_edges,metric_index,c_vec,activator_flag)
    edge_key = [7 8 5 6 3 4 1 2];
    sym_edges = edge_key(~ismember(edge_key,opt_edges));
    rate_vec = NaN(1:8);
    rate_vec(opt_edges) = opt_rates;
    rate_vec(sym_edges) = rate_vec(edge_key(sym_edges));
    metric_vec = calculateMetricsCore(rate_vec,c_vec,activator_flag);
    metric = metric_vec(metric_index);
    