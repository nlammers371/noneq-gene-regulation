function [metric_vec, metric_names, metric_ub_vec, metric_lb_vec] = ...
                              calculateMetricsSym(param_array, simInfo)                            
                            
if ~isempty(param_array)              
    
    % generate input cell     
    paramCell = mat2cell(param_array,size(param_array,1),ones(1,size(param_array,2)));

    % check for saturating behavior
    param_array_high = param_array;
    param_array_high(:,simInfo.cr_index) = 1e2;
    paramCellHigh = mat2cell(param_array_high,size(param_array_high,1),...
                                          ones(1,size(param_array_high,2)));
    
    ProductionRateHigh = productionRateFunction(paramCellHigh{:});

    % calculate production rate         
    ProductionRate = productionRateFunction(paramCell{:});
    
    % calculate sharpness     
    Sharpness = sharpnessFunction(paramCell{:})*simInfo.crs;
    SharpnessNormed = Sharpness./(ProductionRate.*(1-ProductionRate));
    
    % calculate variance    
    Variance = intrinsicVarianceFunction(paramCell{:});
    
    % get state probabilities
    stateProbs = steadyStateVecFunction(paramCell{:});        

    % get state prob entropy
    stateEntropyVec = sum(stateProbs.*log(stateProbs),2);
    stateEntropyVec = stateEntropyVec - min(stateEntropyVec);
    stateEntropyVec = stateEntropyVec/nanmax(stateEntropyVec);

    if simInfo.nStates == 6
        p3Right = stateProbs(:,3);
        p3Wrong = stateProbs(:,5);
    else
        p3Right = NaN(size(ProductionRate));
        p3Wrong = NaN(size(ProductionRate));
    end
    % generate pseudo four state network to calculate specificity
    if simInfo.nStates > 4              
                  
        % extract cw values
        cwVec = param_array(:,simInfo.cw_index);
        
        
        fourFlags = simInfo.fourStateInputFlags;        
        valCellCSFour = mat2cell(param_array(:,fourFlags),size(param_array,1),ones(1,size(param_array(:,fourFlags),2)));    
        
        % right cycle 4 state        
        ProductionRateRight = productionRateFunctionFourState(valCellCSFour{:});
        
        SharpnessRight = sharpnessFunctionFourState(valCellCSFour{:})*simInfo.crs;
        SharpnessRightNorm = SharpnessRight ./ (ProductionRateRight.*(1-ProductionRateRight));
       
        % wrong cycle
        fourWrongFlags = simInfo.fourStateWrongInputFlags;
        valCellCSFourWrong = mat2cell(param_array(:,fourWrongFlags),size(param_array,1),ones(1,size(param_array(:,fourWrongFlags),2)));    
        ProductionRateWrong = productionRateWrongFunctionFourState(valCellCSFourWrong{:});
        
        % inter-network
        specFactorAlt = log10((ProductionRateRight./ProductionRateWrong)./(simInfo.specFactor/simInfo.cw));
        deviationFactor = sharpnessWrongFunction(paramCell{:});
        deviationFactor = 1./deviationFactor;
        
        % intra-network             
        right_flags = simInfo.n_right_bound.*simInfo.activeStateFilter;
        wrong_flags = simInfo.n_wrong_bound.*simInfo.activeStateFilter;
        right_weights = sum(right_flags.*stateProbs,2);
        wrong_weights = sum(wrong_flags.*stateProbs,2);
        specificityFactor = log10((right_weights./wrong_weights)./(simInfo.specFactor./cwVec));                              
        fe_drop = NaN(size(specificityFactor));
        rate_entropy_vec = NaN(size(specificityFactor));
        
    else
        % Free energy drop
        fe_drop = log(prod(param_array(:,simInfo.forwardRateIndices{1}),2)./prod(param_array(:,simInfo.backwardRateIndices{1}),2));
        % Cycle rate J        
%         J = stateProbs(:,1).*param_array(:,simInfo.forwardRateIndices{1}(1)) - ...
%               stateProbs(:,2).*param_array(:,simInfo.backwardRateIndices{1}(1));
        
        
        ProductionRateRight = NaN(size(ProductionRate));
        specificityFactor = NaN(size(ProductionRate));
        specFactorAlt = NaN(size(ProductionRate));
        SharpnessRight = specificityFactor; 
        SharpnessRightNorm = specificityFactor;
        deviationFactor = specificityFactor;
        cwVec = specificityFactor;

        rate_array_norm = param_array(:,2:end)./ sum(param_array(:,2:end),2);
        rate_entropy_vec = sum(rate_array_norm.*log(rate_array_norm),2);
        rate_entropy_vec = rate_entropy_vec - min(rate_entropy_vec);
        rate_entropy_vec = rate_entropy_vec/nanmax(rate_entropy_vec);
    end
    
    % calculate decision metrics 
    [V,T,R] = calculateDecisionMetrics(param_array,simInfo);
    
    % calculate cycle time    
    TauOn = TauONFunction(paramCell{:});
    TauOff = TauOFFFunction(paramCell{:});
    TauCycle = TauOff+TauOn;        
    
    % "Affinity" vec
    AffinityVec = log10(TauOff./TauOn);
    
    % calculate entropy rate
    Phi = entropyRateFunction(paramCell{:}).*TauCycle;    
    
    % simple binomial variance
    BinomialVariance = ProductionRate.*(1-ProductionRate);
    
    % generate vector to output
    if simInfo.TauCycleLimit == 0
        PrecisionOut = log(sqrt(TauCycle./Variance));
      
    elseif simInfo.TauCycleLimit > 0
        PrecisionOut = log(sqrt(simInfo.TauCycleLimit./Variance));
      
    end
      
    metric_vec = [ProductionRate, Sharpness, SharpnessNormed, fe_drop, PrecisionOut,...                
                  Phi, TauCycle,specificityFactor,deviationFactor,ProductionRateRight,...                            
                  SharpnessRight, SharpnessRightNorm, ...
                  log(TauOn./TauOff),V.*TauCycle,1./T,R,...                  
                  log(1./(BinomialVariance)),TauOn,TauOff,specFactorAlt,AffinityVec,log10(cwVec),stateEntropyVec,rate_entropy_vec,p3Right,p3Wrong];
        
    im_flags =  any(round(real(metric_vec),4)~=round(metric_vec,4) & ~isnan(metric_vec),2);
    metric_vec(im_flags,:) = NaN;
    metric_vec = real(metric_vec);
    
    % apply flux limit (default is Inf)
    metric_vec(Phi > simInfo.PhiLimit,:) = NaN;
    metric_vec(TauCycle < simInfo.TauCycleLimit,:) = NaN;
    if simInfo.saturationFlag==1
        metric_vec(round(ProductionRateHigh,3)<0.99,:) = NaN;
    end
else
    metric_vec = [];  
end

metric_names = {'Production Rate','Sharpness','SharpnessNormed','Flux','Precision',...
                'Phi', 'CycleTime','Specificity','deviationFactor','ProductionRateRight',...
                'SharpnessRight','SharpnessRightNorm','ONOFFRatio',...
                'DecisionRateNorm','DecisionTimeNorm','ProductionBinary',...
                'BinomialNoise','TauOn','TauOff','specFactorAlt','AffinityVec','CW','stateEntropy','rateEntropy','p3Right','p3Wrong'};

% specify default bounds
metric_ub_vec = repelem(Inf,length(metric_names));
metric_lb_vec = repelem(-Inf,length(metric_names));
