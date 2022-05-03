function [metric_vec, valArrayCS, metric_names, metric_ub_vec, metric_lb_vec] = ...
                              calculateMetricsNumeric_v3_temp(valArrayCS, simInfo)                            

metric_names = {'Production Rate','Sharpness','Flux','Precision',...
  'Phi', 'CycleTime','Specificity','CWSharpness','deviationFactor',...
  'SharpnessRight', 'SharpnessRightNorm',...
  'DecisionRateNorm','DecisionTimeNorm',...
  'DecisionRate','BinomialNoise','TauOn','TauOff','specFactorAlt','specificityFactorFull','AffinityVec','CW'};

if ~isempty(valArrayCS)      
  
    numerical_precision = simInfo.numerical_precision;        
    
    % initialize parameter arrays
    state_probs_c0 = NaN(size(valArrayCS,1),simInfo.nStates);  
    state_probs_cs = NaN(size(valArrayCS,1),simInfo.nStates);  
    state_probs_c1 = NaN(size(valArrayCS,1),simInfo.nStates);  
    Variance = NaN(size(valArrayCS,1),1);
    VarianceC0 = NaN(size(valArrayCS,1),1);
    VarianceC1 = NaN(size(valArrayCS,1),1);
    TauCycle = NaN(size(valArrayCS,1),1);
    TauOn = NaN(size(valArrayCS,1),1);
    TauOff = NaN(size(valArrayCS,1),1);
    Phi = NaN(size(valArrayCS,1),1);
    
    if simInfo.nStates > 4  && ~simInfo.specOnlyFlag
        % initialize arrays if necessary
        state_probs_c0_right = NaN(size(valArrayCS,1),simInfo.nStates-sum(simInfo.n_wrong_bound>0));  
        state_probs_cs_right = NaN(size(valArrayCS,1),simInfo.nStates-sum(simInfo.n_wrong_bound>0));  
        state_probs_c1_right = NaN(size(valArrayCS,1),simInfo.nStates-sum(simInfo.n_wrong_bound>0));  
        
        state_probs_cs_wrong = NaN(size(valArrayCS,1),simInfo.nStates-sum(simInfo.n_right_bound>0));  
    end
    % generate input arrays
%     valArrayCS = param_array;
        
    % "right" factor
    valArrayC0 = valArrayCS;
    valArrayC0(:,simInfo.cr_index) = simInfo.cr0;    
    
    valArrayC1 = valArrayCS;
    valArrayC1(:,simInfo.cr_index) = simInfo.cr1;     
    
    % helper vec
    activeStateFilter = simInfo.activeStateFilter;    
    use_flags = false(size(valArrayCS,1),1);
    % calculate metric values for each row
    for i = 1:size(valArrayCS,1)
        if all(~isnan(valArrayCS(i,:)) & ~isinf(valArrayCS(i,:)))
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

            % calculate center values for variance, entropy rate, and Tau
            Z_num_cs = calculate_Z_matrix(Q_num_cs,state_probs_cs(i,:),numerical_precision);
            Variance(i) = calculate_var_num(Q_num_cs,state_probs_cs(i,:),activeStateFilter,Z_num_cs,numerical_precision);
            Phi(i) = calculate_entropy_rate_num(Q_num_cs,state_probs_cs(i,:),numerical_precision);   
            if Variance(i)>0 && Phi(i)>=0 && sum(state_probs_cs(i,activeStateFilter))<0.999                
                                   
                [TauOn(i),TauOff(i),TauCycle(i)] = calculate_tau_num(Q_num_cs,state_probs_cs(i,:),activeStateFilter,numerical_precision);
                
                % apply magnitude constraints
%                 if isfield(simInfo,'cw_index')
%                     valArrayCS(i,4:end) = valArrayCS(i,4:end).*TauCycle(i)/simInfo.TauCycleLimit;
%                 else
%                     valArrayCS(i,2:end) = valArrayCS(i,2:end).*TauCycle(i)/simInfo.TauCycleLimit;
%                 end
                oob_flag = false;%any(log10(valArrayCS(i,simInfo.sweepFlags))<simInfo.paramBounds(1,simInfo.sweepFlags)...
%                         |log10(valArrayCS(i,simInfo.sweepFlags))>simInfo.paramBounds(2,simInfo.sweepFlags));  
                if oob_flag
                    valArrayCS(i,:) = NaN;
                else                  
                    use_flags(i) = true;
                    
                    % update variance and Phi
%                     Variance(i) = Variance(i).*simInfo.TauCycleLimit./TauCycle(i);
%                     Phi(i) = Phi(i).*TauCycle(i)./simInfo.TauCycleLimit;                   
                    
                    % update c0 and c1 stuff
%                     Q_num_c0 = Q_num_c0.*TauCycle(i)/simInfo.TauCycleLimit;
%                     Q_num_c1 = Q_num_c1.*TauCycle(i)/simInfo.TauCycleLimit;
                    
                    % reset cycle time
%                     TauOn(i) = TauOn(i).*simInfo.TauCycleLimit./TauCycle(i);
%                     TauOff(i) = TauOff(i).*simInfo.TauCycleLimit./TauCycle(i);
%                     TauCycle(i) = simInfo.TauCycleLimit;
                    
                    % variance for high and low concentration
                    Z_num_c0 = calculate_Z_matrix(Q_num_c0,state_probs_c0(i,:),numerical_precision);
                    Z_num_c1 = calculate_Z_matrix(Q_num_c1,state_probs_c1(i,:),numerical_precision);
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
    
    Sharpness = (ProductionRate1-ProductionRate0)./(simInfo.cr1-simInfo.cr0);
    ProductionRate = sum(state_probs_cs(:,activeStateFilter),2);
    
    % calculate "wrong" sharpness
    CWSharpness = NaN(size(ProductionRate));%(ProductionRate1W-ProductionRate0W)./(cw1-cw0) ./ (pm.*(1-pm));
    
    % calculate decision metrics 
    minError = simInfo.minError;
    [V,T,~] = calculateDecisionMetricsNumeric(ProductionRate0,ProductionRate1,VarianceC0,VarianceC1,minError);    
    VNorm = V.*TauCycle;

    % "Affinity" vec
    AffinityVec = log10(TauOff./TauOn);
        
    % perform specificity-related calculations if appropriate    
    if simInfo.nStates > 4  && ~simInfo.specOnlyFlag
        cwVec = valArrayCS(:,simInfo.cw_index);
        n_right_bound =  simInfo.n_right_bound;
        n_wrong_bound =  simInfo.n_wrong_bound;
        activeStateFilterRight = activeStateFilter(n_wrong_bound==0);
        activeStateFilterWrong = activeStateFilter(n_right_bound==0);
                                
        % right cycle      
        ProductionRateRight = sum(state_probs_cs_right(:,activeStateFilterRight),2);
        ProductionRateRightC1 = sum(state_probs_c1_right(:,activeStateFilterRight),2);
        ProductionRateRightC0 = sum(state_probs_c0_right(:,activeStateFilterRight),2);
        
        SharpnessRight = (ProductionRateRightC1-ProductionRateRightC0)./(simInfo.cr1-simInfo.cr0);
        SharpnessRightNorm = SharpnessRight ./ (ProductionRateRight.*(1-ProductionRateRight));
       
        % wrong cycle
        ProductionRateWrong = sum(state_probs_cs_wrong(:,activeStateFilterWrong),2);
        
        % inter-network
        specFactorParallel = log10((ProductionRateRight./ProductionRateWrong)./(simInfo.specFactor/simInfo.cw));
        deviationFactor = sign((ProductionRateRight-ProductionRate)).*log10(abs(ProductionRate./(ProductionRateRight-ProductionRate)));
        
        % intra-network         
        max_r = max(simInfo.n_right_bound);
        max_w = max(simInfo.n_wrong_bound);
        right_flags = (simInfo.n_right_bound==max_r)&simInfo.activeStateFilter;
        wrong_flags = (simInfo.n_wrong_bound==max_w)&simInfo.activeStateFilter;
        right_weights = sum(right_flags.*state_probs_cs,2);
        wrong_weights = sum(wrong_flags.*state_probs_cs,2);  
        specificityFactorAlt = log10((right_weights./wrong_weights)./simInfo.specFactor.*(cwVec.^max_r)); 
        
        % active states only
        right_weights = sum(simInfo.n_right_bound.*state_probs_cs.*double(simInfo.activeStateFilter),2);
        wrong_weights = sum(simInfo.n_wrong_bound.*state_probs_cs.*double(simInfo.activeStateFilter),2);
        specificityFactor = log10((right_weights./wrong_weights)./simInfo.specFactor.*(cwVec)); 
        
    else
        specFactorParallel = NaN(size(Variance));
        specificityFactor = NaN(size(Variance));
        specificityFactorAlt = NaN(size(Variance));
        SharpnessRightNorm = NaN(size(Variance));
        SharpnessRight = NaN(size(Variance));
        deviationFactor = NaN(size(Variance));
        cwVec = NaN(size(Variance));
    end
    
    % simple binomial noise
    BinomialVariance = ProductionRate.*(1-ProductionRate);
    
    % initialize dummy fields for now fore remaining metrics
    fe_drop = NaN(size(Variance));        
    
    % generate vector to output
    metric_vec = [ProductionRate, Sharpness./BinomialVariance, fe_drop, log(sqrt(BinomialVariance.^2./Variance)),...                
                  Phi.*TauCycle, TauCycle,specificityFactor,CWSharpness,deviationFactor,...                            
                  SharpnessRight, SharpnessRightNorm, ...
                  VNorm,TauCycle./T,V,...                 
                  log(1./(BinomialVariance)),TauOn,TauOff,specFactorParallel,specificityFactorAlt,AffinityVec,log10(cwVec)];
                
        
    % deal with imaginary values
    im_flags =  any(round(real(metric_vec(:,simInfo.edge_metric_indices)),4)~=...
                round(metric_vec(:,simInfo.edge_metric_indices),4) & ~isnan(metric_vec(:,simInfo.edge_metric_indices)),2);
    metric_vec(im_flags,:) = NaN;
    metric_vec(~use_flags) = NaN;
    metric_vec = real(metric_vec);

else
    metric_vec = [];  
end



% specify default bounds
metric_ub_vec = repelem(Inf,length(metric_names));
metric_lb_vec = repelem(-Inf,length(metric_names));
