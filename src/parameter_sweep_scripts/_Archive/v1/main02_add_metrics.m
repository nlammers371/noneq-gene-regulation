clear
close all
% specify paths
addpath('../utilities')
readPath = '../../out/parameter_search/';
metric_names = {'cooperativity', 'information_rate', 'critical_time', 'production_sign', 'production_rate'};
% suffix = 'cooperativity_information_rate';
metric_indices = [1 5];
sim_type = 'edge_sim';
    
tic
for i = 1:2
    equilibrium_flag = i-1;
    if strcmpi(sim_type,'mcmc')
        read_name = ['mcmc_opt_results_eq' num2str(equilibrium_flag) '_' metric_names{metric_indices(1)} '.mat'];
        load([readPath read_name])
        simulation_results = mcmc_results;
    else
        read_name = ['edge_sim_results_eq' num2str(equilibrium_flag) '_' ...
            metric_names{metric_indices(1)} '_' metric_names{metric_indices(2)} '.mat'];
        load([readPath read_name])
    end    
    % check to see if energy dissapation info is present
    if ~isfield(simulation_results,'energy_per_cycle')
        for k = 1:numel(simulation_results)
            % calculate dissipation per cycle        
            k_mat = simulation_results(k).rate_array(:,1:4);
            r_mat = simulation_results(k).rate_array(:,5:8);
            e_vec = log(prod(k_mat,2) ./ prod(r_mat,2));
            simulation_results(k).energy_per_cycle = e_vec;
            % calculate net cycle flux
            ss_mat = NaN(numel(e_vec),4);
            flux_vec = NaN(size(e_vec));
            real_flag_vec = false(size(e_vec));
            for j = 1:size(r_mat,1)
                r_vec = r_mat(j,:);
                if any(isnan(r_vec))
                    continue
                end
                k_vec = k_mat(j,:);
                R = [-k_vec(1)-r_vec(4)   r_vec(1)             0               k_vec(4); 
                    k_vec(1)        -r_vec(1)-k_vec(2)      r_vec(2)              0
                    0           k_vec(2)         -r_vec(2)-k_vec(3)         r_vec(3)
                    r_vec(4)          0               k_vec(3)        -r_vec(3)-k_vec(4) ];
                [V,D] = eig(R);
                [~,mi] = max(diag(real(D)));
                ss_vec = V(:,mi)/sum(V(:,mi));
                ss_mat(j,:) = ss_vec;
                flux =  ss_vec(1)*k_vec(1) - ss_vec(2)*r_vec(1);     
                real_flag_vec(j) = isreal(D);               
                flux_vec(j) = flux;
            end
            simulation_results(k).ss_mat = ss_mat;
            simulation_results(k).real_flag_vec = real_flag_vec;
            simulation_results(k).cycle_flux_vec = flux_vec;
        end
    end
    if strcmpi(sim_type,'mcmc')
        mcmc_results = simulation_results;
        save([readPath read_name],'mcmc_results')
    else
        save([readPath read_name],'simulation_results')
    end
end
toc