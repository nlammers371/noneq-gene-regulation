function s0_f0_motif_struct = generate_s_vs_f0_vs_mu_scatters(nBins,sweep_options,beta,FigPath)

    [~,metric_names] = calculateMetricsMultiState([]);
    % get sweep indices
    spec_index = find(strcmp(metric_names,'Specificity'));    
    sharp_right_index = find(strcmp(metric_names,'SharpnessRight'));
    sharp_right_norm_index = find(strcmp(metric_names,'SharpnessRightNorm'));    
    sharpness_index = find(strcmp(metric_names,'Sharpness'));
    
    % initialize
    mu_vec = logspace(log10(1),log10(beta^2),nBins);

    s0_f0_motif_struct = struct;
    s0_f0_motif_struct.s_vec_sharp = NaN(1,length(mu_vec));
    s0_f0_motif_struct.s_vec_spec = NaN(1,length(mu_vec));
    s0_f0_motif_struct.spec_vec = NaN(1,length(mu_vec));
    s0_f0_motif_struct.s0_vec = NaN(1,length(mu_vec));
    s0_f0_motif_struct.s_max_vec = NaN(1,length(mu_vec));
    s0_f0_motif_struct.mu_vec = mu_vec;
    
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
        [~, sim_struct_sharp] = param_sweep_multi([sharp_right_norm_index sharpness_index],sweep_options{:},...
                                                  'half_max_flag',true,'wrongFactorConcentration',mu_vec(m),'equilibrium_flag',...
                                                  false,'n_sim',5,'useParpool',0,'specFactor',beta);

        [~, sim_struct_spec] = param_sweep_multi([spec_index sharpness_index],sweep_options{:},...
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
        
        [s0_f0_motif_struct(m).s_vec_sharp, mi_sharp] = nanmax(result_array_sharp(sharp_ind,sharpness_index));    
        s0_f0_motif_struct(m).s0_vec = result_array_sharp(sharp_ind(mi_sharp),sharp_right_norm_index);         
        [s0_f0_motif_struct(m).s_vec_spec, mi_spec] = nanmax(result_array_spec(spec_ind,sharpness_index));    
        s0_f0_motif_struct(m).spec_vec = result_array_spec(spec_ind(mi_spec),spec_index);

        s0_f0_motif_struct(m).s_max_vec = nanmax(vertcat(result_array_spec(:,sharpness_index),result_array_sharp(:,sharpness_index)));
        toc   
    end                                            

    save([FigPath 's0_f0_motif_struct.mat'],'s0_f0_motif_struct')