function [metric_vec, param_array_out, metric_names, metric_ub_vec, metric_lb_vec] = ...
                              calculateMetricsNumeric_v3(param_array, simInfo)                            

metric_names = {'ProductionRate','Sharpness','Precision','Phi','IR',...
                'Specificity','SharpnessRight','TauCycle',...                
                'SharpnessRaw','VarianceRaw','SharpnessRightRaw',...
                'InverseDecisionTime','TauOn','TauOff','specFactorParallel','CW'};
              
% specify default bounds   
metric_ub_vec = repelem(Inf,length(metric_names));
metric_lb_vec = repelem(-Inf,length(metric_names));

% initialize output param array
param_array_out = param_array;                            
if ~isempty(param_array)      
  
    % set matrix near singular warning to error status
    warning('error', 'MATLAB:nearlySingularMatrix');
    
    % NL: these bounds prevent common errors related to numerical precision
    phi_index = find(strcmp(metric_names,'Phi'));
    cw_index = find(strcmp(metric_names,'CW'));
    tau_index = strcmp(metric_names,'TauCycle'); 
    precision_index = strcmp(metric_names,'Precision'); 
    sharpness_index = strcmp(metric_names,'Sharpness'); 
    if ~simInfo.equilibrium_flag
        if any(ismember(simInfo.edge_metric_indices,[phi_index]))
            metric_lb_vec(phi_index) = 1e-3;        
        else
            metric_lb_vec(phi_index) = 1e-2; % NL: this seems to avoid most precision errors
        end
        if any(ismember(simInfo.edge_metric_indices,[cw_index]))
            metric_ub_vec(phi_index) = 5e5;        
        else
            metric_ub_vec(phi_index) = 5e4; 
        end        
    end    
    metric_lb_vec(sharpness_index) = 1e-3;    
    metric_lb_vec(precision_index) = 1e-3;
    metric_ub_vec(tau_index) = 5e4;

    % load other relevant parameters into working space
    numerical_precision = simInfo.numerical_precision;        
    w_flags = contains(simInfo.sweepVarStrings,'w')|strcmp(simInfo.sweepVarStrings,'a')|strcmp(simInfo.sweepVarStrings,'cw');
    sweepFlags = simInfo.sweepFlags;  
    paramBounds = simInfo.paramBounds;  
    paramBoundsRel = paramBounds(:,sweepFlags&~w_flags); % indicates variables to adjust with tau_c
    
    % extract upper and lower bounds to use for expression rate
    a1 = simInfo.a1;
    a0 = simInfo.a0;
    
    state_probs_c0 = NaN(size(param_array,1),simInfo.nStates);  
    state_probs_cs = NaN(size(param_array,1),simInfo.nStates);  
    state_probs_c1 = NaN(size(param_array,1),simInfo.nStates);  
    VarianceRaw = NaN(size(param_array,1),1);
    VarianceC0 = NaN(size(param_array,1),1);
    VarianceC1 = NaN(size(param_array,1),1);
    TauCycle = NaN(size(param_array,1),1);
    TauOn = NaN(size(param_array,1),1);
    TauOff = NaN(size(param_array,1),1);
    Phi = NaN(size(param_array,1),1);
    
    if simInfo.nStates > 4  && ~simInfo.specOnlyFlag
        % initialize arrays if necessary
        state_probs_c0_right = NaN(size(param_array,1),simInfo.nStates-sum(simInfo.n_wrong_bound>0));  
        state_probs_cs_right = NaN(size(param_array,1),simInfo.nStates-sum(simInfo.n_wrong_bound>0));  
        state_probs_c1_right = NaN(size(param_array,1),simInfo.nStates-sum(simInfo.n_wrong_bound>0));          
        state_probs_cs_wrong = NaN(size(param_array,1),simInfo.nStates-sum(simInfo.n_right_bound>0));  
    end    
        
    % initialize arrays    
    valArrayC0 = NaN(size(param_array));
    valArrayC1 = NaN(size(param_array));
    valArrayCS = NaN(size(param_array));
    
    % helper vec that indicates active states
    activeStateFilter = simInfo.activeStateFilter;
    
    use_flags = false(size(param_array,1),1);
    % calculate metric values for each row
    for i = 1:size(param_array,1)
        if all(~isnan(param_array(i,:)) & ~isinf(param_array(i,:)))        

            valCellInit = mat2cell(param_array(i,:),size(param_array(i,:),1),ones(1,size(param_array(i,:),2)));    

            % get rate arrays
            Q_num_cs = RSymFun(valCellInit{:});            
            cs_probs_temp = calculate_ss_num(Q_num_cs,numerical_precision);    

            % calculate cycle time (will use for normalization)
            try
                [TauOn(i),TauOff(i),TauCycle(i)] = ...
                        calculate_tau_num(Q_num_cs,cs_probs_temp,activeStateFilter,numerical_precision);
            catch
                % do nothing               
            end
            
            % normalize rate values and check that they are in bounds
            params_new = param_array(i,:);
            params_new(sweepFlags&~w_flags) = params_new(sweepFlags&~w_flags)*TauCycle(i)./simInfo.TauCycleTime;

            % check to see if anything has been pushed out of range
            oob_flag = max(log10(params_new(sweepFlags&~w_flags))<paramBoundsRel(1,:)...
                        |log10(params_new(sweepFlags&~w_flags))>paramBoundsRel(2,:),[],2);  
                   
            oob_flag = oob_flag | any(isnan(params_new));

            if ~oob_flag
                valArrayCS(i,:) = params_new;
                valArrayC0(i,:) = params_new;
                valArrayC0(i,simInfo.cr_index) = simInfo.cr0;
                valArrayC1(i,:) = params_new;
                valArrayC1(i,simInfo.cr_index) = simInfo.cr1;

                valCellCS = mat2cell(valArrayCS(i,:),size(valArrayCS(i,:),1),ones(1,size(valArrayCS(i,:),2)));    
                valCellC0 = mat2cell(valArrayC0(i,:),size(valArrayC0(i,:),1),ones(1,size(valArrayC0(i,:),2)));
                valCellC1 = mat2cell(valArrayC1(i,:),size(valArrayC1(i,:),1),ones(1,size(valArrayC1(i,:),2)));

                % get rate arrays
                Q_num_cs = RSymFun(valCellCS{:});
                Q_num_c0 = RSymFun(valCellC0{:});
                Q_num_c1 = RSymFun(valCellC1{:});

                % calculate probabilities                
                state_probs_cs(i,:) = calculate_ss_num(Q_num_cs,numerical_precision);                      
                state_probs_c0(i,:) = calculate_ss_num(Q_num_c0,numerical_precision);
                state_probs_c1(i,:) = calculate_ss_num(Q_num_c1,numerical_precision);
                
%                 low_flag = false;%log10(nanmin(state_probs_cs(i,:))) <= -numerical_precision;
%                 if all(ismember(simInfo.edge_metric_indices,[precision_index, sharpness_index]))
%                     low_flag = false;%log10(nanmin(state_probs_cs(i,:))) <= -numerical_precision;
%                 end
                % calculate center values for variance, entropy rate, and Tau
                try % NL: this calculation requres a matrix inversion. Will lead to error if matrix is near singular
                    Z_num_cs = calculate_Z_matrix(Q_num_cs,state_probs_cs(i,:),numerical_precision);
                catch
                    Z_num_cs = NaN(size(Q_num_cs));
                end
                VarianceRaw(i) = calculate_var_num(Q_num_cs,state_probs_cs(i,:),activeStateFilter,Z_num_cs,numerical_precision);
                Phi(i) = calculate_entropy_rate_num(Q_num_cs,state_probs_cs(i,:),numerical_precision);   

                if ~low_flag && VarianceRaw(i)>0 && Phi(i)>=0 && sum(state_probs_cs(i,activeStateFilter))<a1 && ...
                                                  sum(state_probs_cs(i,activeStateFilter))>a0  
                    use_flags(i) = true;

                    % variance for high and low concentration
                    try
                        Z_num_c0 = calculate_Z_matrix(Q_num_c0,state_probs_c0(i,:),numerical_precision);
                        Z_num_c1 = calculate_Z_matrix(Q_num_c1,state_probs_c1(i,:),numerical_precision);
                    catch
                        Z_num_c0 = NaN(size(Q_num_c0));
                        Z_num_c1 = NaN(size(Q_num_c1));
                    end
                    VarianceC0(i) = calculate_var_num(Q_num_c0,state_probs_c0(i,:),activeStateFilter,Z_num_c0,numerical_precision);
                    VarianceC1(i) = calculate_var_num(Q_num_c1,state_probs_c1(i,:),activeStateFilter,Z_num_c1,numerical_precision);

                    % perform specificity-related calculations if appropriate
                    if simInfo.nStates > 4  && ~simInfo.specOnlyFlag                
                        n_wrong_bound =  simInfo.n_wrong_bound;
                        n_right_bound =  simInfo.n_right_bound;

                        % truncate numerical arrays
                        Q_num_cs_right = Q_num_cs(n_wrong_bound==0,n_wrong_bound==0);
                        Q_num_cs_right(eye(size(Q_num_cs_right,1))==1) = 0;
                        Q_num_cs_right(eye(size(Q_num_cs_right,1))==1) = -sum(Q_num_cs_right);

                        Q_num_c1_right = Q_num_c1(n_wrong_bound==0,n_wrong_bound==0);
                        Q_num_c1_right(eye(size(Q_num_c1_right,1))==1) = 0;
                        Q_num_c1_right(eye(size(Q_num_c1_right,1))==1) = -sum(Q_num_c1_right);

                        Q_num_c0_right = Q_num_c0(n_wrong_bound==0,n_wrong_bound==0);
                        Q_num_c0_right(eye(size(Q_num_c0_right,1))==1) = 0;
                        Q_num_c0_right(eye(size(Q_num_c0_right,1))==1) = -sum(Q_num_c0_right);

                        % perform basic calculations using "right" network
                        % calculate probabilities
                        state_probs_cs_right(i,:) = calculate_ss_num(Q_num_cs_right,numerical_precision);        
                        state_probs_c0_right(i,:) = calculate_ss_num(Q_num_c0_right,numerical_precision);
                        state_probs_c1_right(i,:) = calculate_ss_num(Q_num_c1_right,numerical_precision);

                        % generate equivalent "wrong" network
                        Q_num_cs_wrong = Q_num_cs(n_right_bound==0,n_right_bound==0);
                        Q_num_cs_wrong(eye(size(Q_num_cs_wrong,1))==1) = 0;
                        Q_num_cs_wrong(eye(size(Q_num_cs_wrong,1))==1) = -sum(Q_num_cs_wrong);
                        state_probs_cs_wrong(i,:) = calculate_ss_num(Q_num_cs_wrong,numerical_precision);                 
                    end
                end
            end
        end
    end    

    % calculate sharpness and production rate
    ProductionRate1 = sum(state_probs_c1(:,activeStateFilter),2);
    ProductionRate0 = sum(state_probs_c0(:,activeStateFilter),2);   
    ProductionRate = sum(state_probs_cs(:,activeStateFilter),2);
    
    % calculate sharpness 3 different ways and compare to check for
    % numerical precision errors in SS vectors    
    dp10 = ProductionRate1-ProductionRate0;
    SharpnessRaw01 = dp10./(simInfo.cr1-simInfo.cr0);   
    SharpnessRawS0 = (ProductionRate-ProductionRate0)./(simInfo.crs-simInfo.cr0);   
    SharpnessRaw1S = (ProductionRate1-ProductionRate)./(simInfo.cr1-simInfo.crs);   
    
    % calculate b terms
    b_factor = (ProductionRate.*(1-ProductionRate));
    pd1S = (ProductionRate1+ProductionRate)/2;
    b_factor1S = (pd1S.*(1-pd1S));
    pdS0 = (ProductionRate0+ProductionRate)/2;
    b_factorS0 = (pdS0.*(1-pdS0));
    
    % calculate normalized sharpness terms
    Sharpness = SharpnessRaw01 ./ b_factor .* simInfo.crs;
    SharpnessS0 = SharpnessRawS0 ./ b_factorS0 .* (simInfo.crs+simInfo.cr0)/2;
    Sharpness1S = SharpnessRaw1S ./ b_factor1S .* (simInfo.crs+simInfo.cr1)/2;
    
    % look for unusual deviations
    s_array = [SharpnessS0 Sharpness1S Sharpness];
    sig = nanstd(real(s_array),[],2); % note that we screen separately for imaginary components
    mu = nanmin(real(s_array),[],2);
    
    % NL: 0.1 threshold is ad hoc, but seems to do a good job
    % picking out erroneous cases arising from precision issues
    qc_flags = abs(sig./mu) <= 0.1 & ~any(isnan(s_array),2) & all(real(s_array)>=0,2);%all(round(real(Sharpness),1) == round(real([Sharpness SharpnessS0 Sharpness1S]),1),2);
    if all(ismember(simInfo.edge_metric_indices,[precision_index, sharpness_index]))
        use_flags = use_flags & qc_flags; % designed specfically to deal with artifactis in S vs P sweeps
    end
    
%     qc_flag_alt = round(SharpnessRaw01,4) == round(mean(abs([SharpnessRawS0 SharpnessRaw1S]),2),4)&all(sign(SharpnessRaw01)==sign([SharpnessRawS0 SharpnessRaw1S]),2);
    Variance = VarianceRaw ./ b_factor.^2;
    VarianceC0Norm = VarianceC0 ./ (ProductionRate0.*(1-ProductionRate0)).^2;
    VarianceC1Norm = VarianceC1 ./ (ProductionRate1.*(1-ProductionRate1)).^2;
    
    v_array = [Variance VarianceC1Norm VarianceC0Norm];
    sig_var = nanstd(real(v_array),[],2); % note that we screen separately for imaginary components
    mu_var = nanmin(real(v_array),[],2);
    qc_flags_var = sig_var./mu_var <= .4 & all(v_array>0,2);
    if all(ismember(simInfo.edge_metric_indices,[precision_index, sharpness_index]))
        use_flags = use_flags & qc_flags_var;
    end
    % calculate "wrong" sharpness
%     CWSharpness = NaN(size(ProductionRate));%(ProductionRate1W-ProductionRate0W)./(cw1-cw0) ./ (pm.*(1-pm));
    
    % calculate decision metrics 
    minError = simInfo.minError;
    [V,T,~] = calculateDecisionMetricsNumeric(ProductionRate0,ProductionRate1,VarianceC0,VarianceC1,minError);        

    % perform specificity-related calculations if appropriate    
    if simInfo.nStates > 4  && ~simInfo.specOnlyFlag
        cwVec = param_array(:,simInfo.cw_index);
        n_right_bound =  simInfo.n_right_bound;
        n_wrong_bound =  simInfo.n_wrong_bound;
        activeStateFilterRight = activeStateFilter(n_wrong_bound==0);
        activeStateFilterWrong = activeStateFilter(n_right_bound==0);
                                
        % right cycle      
        ProductionRateRight = sum(state_probs_cs_right(:,activeStateFilterRight),2);
        ProductionRateRightC1 = sum(state_probs_c1_right(:,activeStateFilterRight),2);
        ProductionRateRightC0 = sum(state_probs_c0_right(:,activeStateFilterRight),2);
        
        SharpnessRightRaw = (ProductionRateRightC1-ProductionRateRightC0)./(simInfo.cr1-simInfo.cr0);
        SharpnessRight = SharpnessRightRaw ./ (ProductionRateRight.*(1-ProductionRateRight)) .* simInfo.crs;
       
        % wrong cycle
        ProductionRateWrong = sum(state_probs_cs_wrong(:,activeStateFilterWrong),2);
        
        % inter-network
        specFactorParallel = log10((ProductionRateRight./ProductionRateWrong)./(simInfo.specFactor/simInfo.cw));
        
        % intra-network         
%         max_r = max(simInfo.n_right_bound);
%         max_w = max(simInfo.n_wrong_bound);
%         right_flags = (simInfo.n_right_bound==max_r)&simInfo.activeStateFilter;
%         wrong_flags = (simInfo.n_wrong_bound==max_w)&simInfo.activeStateFilter;
%         right_weights = sum(simInfo.n_right_bound.*state_probs_cs,2);
%         wrong_weights = sum(simInfo.n_wrong_bound.*state_probs_cs,2);  
%         specificityFactorAlt = log10((right_weights./wrong_weights)./simInfo.specFactor.*(cwVec.^max_r)); 
        
        % active states only
        right_weights = sum(simInfo.n_right_bound.*state_probs_cs.*double(simInfo.activeStateFilter),2);
        wrong_weights = sum(simInfo.n_wrong_bound.*state_probs_cs.*double(simInfo.activeStateFilter),2);
        specificityFactor = log10((right_weights./wrong_weights)./simInfo.specFactor.*(cwVec)); 
        
    else
        specFactorParallel = NaN(size(VarianceRaw));
        specificityFactor = NaN(size(VarianceRaw));        
        SharpnessRight = NaN(size(VarianceRaw));
        SharpnessRightRaw = NaN(size(VarianceRaw));        
        cwVec = NaN(size(VarianceRaw));
    end    
    
    % generate vector to output
    metric_vec = [ProductionRate, Sharpness, sqrt(1./Variance),Phi,V,...                
                   specificityFactor,SharpnessRight,TauCycle,...                        
                   SharpnessRaw01,VarianceRaw,SharpnessRightRaw, ...
                   1./T,TauOn,TauOff,specFactorParallel,log10(cwVec)];

        
    % deal with imaginary values
    im_flags =  any(round(real(metric_vec(:,simInfo.edge_metric_indices)),4)~=...
                round(metric_vec(:,simInfo.edge_metric_indices),4) & ~isnan(metric_vec(:,simInfo.edge_metric_indices)),2);
    metric_vec(im_flags,:) = NaN;
    metric_vec(~use_flags,:) = NaN;
    metric_vec = real(metric_vec);
    
    % check for metrics that are oob
    m_oob_flags = any(metric_vec>metric_ub_vec,2) | any(metric_vec<metric_lb_vec,2);
    metric_vec(m_oob_flags,:) = NaN;
    
    % update parameter array
    param_array_out(use_flags,:) = valArrayCS(use_flags,:);
    param_array_out(~use_flags|m_oob_flags,:) = NaN;

else
    metric_vec = [];  
end

