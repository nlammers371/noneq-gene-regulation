function [metric_vec, metric_names, metric_ub_vec, metric_lb_vec] = ...
                              calculateMetricsNumeric_v2(param_array, simInfo)                            
                            
if ~isempty(param_array)      
  
    numerical_precision = simInfo.numerical_precision;        
    
    % initialize parameter arrays
    state_probs_c0 = NaN(size(param_array,1),simInfo.nStates);  
    state_probs_cs = NaN(size(param_array,1),simInfo.nStates);  
    state_probs_c1 = NaN(size(param_array,1),simInfo.nStates);  
    Variance = NaN(size(param_array,1),1);
    VarianceC0 = NaN(size(param_array,1),1);
    VarianceC1 = NaN(size(param_array,1),1);
    TauCycle = NaN(size(param_array,1),1);
    TauOn = NaN(size(param_array,1),1);
    TauOff = NaN(size(param_array,1),1);
    Phi = NaN(size(param_array,1),1);
    
    % generate input arrays
    valArrayCS = param_array;
        
    % "right" factor
    valArrayC0 = valArrayCS;
    valArrayC0(:,simInfo.cr_index) = simInfo.cr0;    
    
    valArrayC1 = valArrayCS;
    valArrayC1(:,simInfo.cr_index) = simInfo.cr1;     
    
    % helper vec
    activeStateFilter = simInfo.activeStateFilter;
    
    % initialize cell array to store
    
    % calculate metric values for each row
    for i = 1:size(param_array,1)
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
            [TauOn(i),TauOff(i),TauCycle(i)] = calculate_tau_num(Q_num_cs,state_probs_cs(i,:),activeStateFilter,numerical_precision);

            % variance for high and low concentration
            Z_num_c0 = calculate_Z_matrix(Q_num_c0,state_probs_c0(i,:),numerical_precision);
            Z_num_c1 = calculate_Z_matrix(Q_num_c1,state_probs_c1(i,:),numerical_precision);
            VarianceC0(i) = calculate_var_num(Q_num_c0,state_probs_c0(i,:),activeStateFilter,Z_num_c0,numerical_precision);
            VarianceC1(i) = calculate_var_num(Q_num_c1,state_probs_c1(i,:),activeStateFilter,Z_num_c1,numerical_precision);
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
        
    % perform calculations for "right cycle" 
    if simInfo.nStates > 4  && ~simInfo.specOnlyFlag
        % assign CW vec values
        cwVec = param_array(:,simInfo.cw_index);
        
        % calculate "locus-centric" specificity
        max_r = max(simInfo.n_right_bound);
        max_w = max(simInfo.n_wrong_bound);
        right_flags = (simInfo.n_right_bound==max_r)&simInfo.activeStateFilter;
        wrong_flags = (simInfo.n_wrong_bound==max_w)&simInfo.activeStateFilter;
        right_weights = sum(right_flags.*state_probs_cs,2);
        wrong_weights = sum(wrong_flags.*state_probs_cs,2);
        specificityFactor = log10((right_weights./wrong_weights)./(simInfo.specFactor./cwVec));
        
        % calculate metrics for "right" and "wrong" 4 state models,
        % mimicking "paralell locus" formulation of specificity
        fourFlags = simInfo.fourStateInputFlags;
        valCellCSFour = mat2cell(valArrayCS(:,fourFlags),size(valArrayCS,1),ones(1,size(valArrayCS(:,fourFlags),2)));    
        valCellC0Four = mat2cell(valArrayC0(:,fourFlags),size(valArrayC0,1),ones(1,size(valArrayC0(:,fourFlags),2)));
        valCellC1Four = mat2cell(valArrayC1(:,fourFlags),size(valArrayC1,1),ones(1,size(valArrayC1(:,fourFlags),2)));
        
        % right cycle 4 state        
        ProductionRateRight = productionRateFunctionFourState(valCellCSFour{:});
        ProductionRateRightC1 = productionRateFunctionFourState(valCellC1Four{:});
        ProductionRateRightC0 = productionRateFunctionFourState(valCellC0Four{:});
        
        SharpnessRight = (ProductionRateRightC1-ProductionRateRightC0)./(simInfo.cr1-simInfo.cr0);
        SharpnessRightNorm = SharpnessRight ./ (ProductionRateRight.*(1-ProductionRateRight));
       
        % wrong cycle
        fourWrongFlags = simInfo.fourStateWrongInputFlags;        
        valCellCSFourWrong = mat2cell(valArrayCS(:,fourWrongFlags),size(valArrayCS,1),ones(1,size(valArrayCS(:,fourWrongFlags),2)));    
        ProductionRateWrong = productionRateWrongFunctionFourState(valCellCSFourWrong{:});
        
        % inter-network
        specFactorAlt = log10((ProductionRateRight./ProductionRateWrong)./(simInfo.specFactor/simInfo.cw));
        deviationFactor = sign((ProductionRateRight-ProductionRate)).*log10(abs(ProductionRate./(ProductionRateRight-ProductionRate)));                 
    else
        specFactorAlt = NaN(size(Variance));
        specificityFactor = NaN(size(Variance));
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
    metric_vec = [ProductionRate, Sharpness, fe_drop, log(sqrt(TauCycle./Variance)),...                
                  Phi.*TauCycle, TauCycle,specificityFactor,CWSharpness,deviationFactor,...                            
                  SharpnessRight, SharpnessRightNorm, ...
                  VNorm,TauCycle./T,V,...                 
                  log(1./(BinomialVariance)),TauOn,TauOff,specFactorAlt,AffinityVec,log10(cwVec)];
                
    % deal with imaginary values
    im_flags =  any(round(real(metric_vec(:,simInfo.edge_metric_indices)),4)~=...
                round(metric_vec(:,simInfo.edge_metric_indices),4) & ~isnan(metric_vec(:,simInfo.edge_metric_indices)),2);
    metric_vec(im_flags,:) = NaN;
    metric_vec = real(metric_vec);

else
    metric_vec = [];  
end

metric_names = {'Production Rate','Sharpness','Flux','Precision',...
                'Phi', 'CycleTime','Specificity','CWSharpness','deviationFactor',...
                'SharpnessRight', 'SharpnessRightNorm',...
                'DecisionRateNorm','DecisionTimeNorm',...
                'DecisionRate','BinomialNoise','TauOn','TauOff','specFactorAlt','AffinityVec','CW'};

% specify default bounds
metric_ub_vec = repelem(Inf,length(metric_names));
metric_lb_vec = repelem(-Inf,length(metric_names));
