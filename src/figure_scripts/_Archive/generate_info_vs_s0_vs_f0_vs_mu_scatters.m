function info_s0_f0_motif_struct = generate_info_vs_s0_vs_f0_vs_mu_scatters(nBins,sweep_options,beta,FigPath)

    [~,metric_names] = calculateMetricsMultiState([]);
    % get sweep indices
    spec_index = find(strcmp(metric_names,'Specificity'));    
    info_index = find(strcmp(metric_names,'DecisionRateNorm'));
    sharp_right_norm_index = find(strcmp(metric_names,'SharpnessRightNorm'));    
%     sharpness_index = find(strcmp(metric_names,'Sharpness'));
    
    % initialize
    mu_vec = logspace(log10(1),log10(beta^2),nBins);

    info_s0_f0_motif_struct = struct;
    
    p = gcp('nocreate');
    myCluster = parcluster('local');
    
    NumWorkers = min([myCluster.NumWorkers 24]);    
    if isempty(p)      
        parpool(NumWorkers);      
    elseif p.NumWorkers~=NumWorkers
      delete(p);      
      parpool(NumWorkers);      
    end
    
    parfor m = 1:length(mu_vec)
        tic
        [~, sim_struct_sharp] = param_sweep_multi([sharp_right_norm_index info_index],sweep_options{:},...
                                                  'half_max_flag',false,'wrongFactorConcentration',mu_vec(m),'equilibrium_flag',...
                                                  false,'n_sim',5,'useParpool',0,'specFactor',beta);

        [~, sim_struct_spec] = param_sweep_multi([spec_index info_index],sweep_options{:},...
                                                  'half_max_flag',true,'wrongFactorConcentration',mu_vec(m),...
                                                  'equilibrium_flag',false,'n_sim',5,'specFactor',beta);

        result_array_sharp = vertcat(sim_struct_sharp.metric_array);
        result_array_spec = vertcat(sim_struct_spec.metric_array);

        sharp_95 = prctile(result_array_sharp(:,sharp_right_norm_index),95);
        sharp_ind = find(result_array_sharp(:,sharp_right_norm_index)>=sharp_95);    
        if isempty(sharp_ind)
          error('uh oh')
        end

        spec_ind = find(result_array_spec(:,spec_index)>=1.9);    
        if isempty(spec_ind)
          error('uh oh')
        end
        
        [info_s0_f0_motif_struct(m).info_vec_sharp, mi_sharp] = nanmax(result_array_sharp(sharp_ind,info_index));    
        info_s0_f0_motif_struct(m).s0_vec = result_array_sharp(sharp_ind(mi_sharp),sharp_right_norm_index);         
        [info_s0_f0_motif_struct(m).info_vec_spec, mi_spec] = nanmax(result_array_spec(spec_ind,info_index));    
        info_s0_f0_motif_struct(m).spec_vec = result_array_spec(spec_ind(mi_spec),spec_index);

        info_s0_f0_motif_struct(m).info_max_vec = nanmax(vertcat(result_array_spec(:,info_index),result_array_sharp(:,info_index)));
        info_s0_f0_motif_struct(m).mu_vec = mu_vec(m);
        toc   
    end                                            

    save([FigPath 'info_s0_f0_motif_struct.mat'],'info_s0_f0_motif_struct')