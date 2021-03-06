function [metric_vec, param_array, metric_names, metric_ub_vec, metric_lb_vec] = ...
                              calculateMetricsSym_v2(param_array, simInfo)                            
                            
if ~isempty(param_array)              
    
    % generate input cell     
    paramCell = mat2cell(param_array,size(param_array,1),ones(1,size(param_array,2)));

  
    w_flags = contains(simInfo.sweepVarStrings,'w')|strcmp(simInfo.sweepVarStrings,'a')|strcmp(simInfo.sweepVarStrings,'cw');
    sweepFlags = simInfo.sweepFlags;  
%     paramBounds = simInfo.paramBounds;  
%     paramBoundsRel = paramBounds(:,sweepFlags&~w_flags); % indicates variables to adjust with tau_c
    
    % calculate cycle time    
    TauOn = TauONFunction(paramCell{:});
    TauOff = TauOFFFunction(paramCell{:});
    TauCycle = TauOff+TauOn; 
        
    % normalize rate values and check that they are in bounds
%     params_new = param_array;
    param_array(:,sweepFlags&~w_flags) = param_array(:,sweepFlags&~w_flags).*TauCycle./simInfo.TauCycleTime;
    
    % renormalize by cycle time
%     warning('Need to adjust this so Tau shift is only applied to base rates')
%     if isfield(simInfo,'cw_index')
%         param_array(:,4:end) = param_array(:,4:end).*TauCycle/simInfo.TauCycleTime;
%     else
%         param_array(:,2:end) = param_array(:,2:end).*TauCycle/simInfo.TauCycleTime;
%     end
    oob_flags = max(log10(param_array(:,simInfo.sweepFlags))<simInfo.paramBounds(1,simInfo.sweepFlags)...
            |log10(param_array(:,simInfo.sweepFlags))>simInfo.paramBounds(2,simInfo.sweepFlags),[],2);  

    param_array(oob_flags,:) = NaN;

    paramCell = mat2cell(param_array,size(param_array,1),ones(1,size(param_array,2)));

    TauOn = TauONFunction(paramCell{:});
    TauOff = TauOFFFunction(paramCell{:});
    TauCycle = TauOff+TauOn;   
        
%     ProductionRateHigh = productionRateFunction(paramCellHigh{:});

    % calculate production rate         
    ProductionRate = productionRateFunction(paramCell{:});
    
    % calculate sharpness     
    SharpnessRaw = sharpnessFunction(paramCell{:})*simInfo.crs;
    Sharpness = SharpnessRaw./(ProductionRate.*(1-ProductionRate));
    
    % calculate variance    
    Variance = intrinsicVarianceFunction(paramCell{:});            
%     Variance = VarianceRaw;            
    
    % get state probabilities
    stateProbs = steadyStateVecFunction(paramCell{:});        

    % get state prob entropy
    stateEntropyVec = sum(stateProbs.*log(stateProbs),2);
    stateEntropyVec = stateEntropyVec - min(stateEntropyVec);
    stateEntropyVec = stateEntropyVec/nanmax(stateEntropyVec);
    
    rate_entropy_array = NaN(size(stateEntropyVec));
    if any(strcmp(simInfo.metric_names(simInfo.edge_metric_indices),'rateEntropy')) || simInfo.calculateRateEntropy
        for p = 1:size(param_array,1)
            paramCellTemp = mat2cell(param_array(p,:),size(param_array(p,:),1),ones(1,size(param_array(p,:),2)));
            RTemp = RSymFun(paramCellTemp{:});
            rate_vec = RTemp(RTemp>0);
            
            rate_vec_norm = rate_vec./ sum(rate_vec);
            rate_entropy_array(p,:) = sum(rate_vec_norm.*log(rate_vec_norm));            
        end
        rate_entropy_array = rate_entropy_array - min(rate_entropy_array);
        rate_entropy_array = rate_entropy_array/nanmax(rate_entropy_array);
        
    end
        
    if simInfo.nStates == 6
        p3Right = stateProbs(:,3);
        p3Wrong = stateProbs(:,5);
    else
        p3Right = NaN(size(ProductionRate));
        p3Wrong = NaN(size(ProductionRate));
    end
    
    rand_vec = rand(size(ProductionRate));
    % generate pseudo four state network to calculate specificity
    if simInfo.nStates > 4  && ~simInfo.specOnlyFlag             
                  
        % extract cw values
        cwVec = param_array(:,simInfo.cw_index);        
        
%         fourFlags = simInfo.fourStateInputFlags;        
%         temp_array = param_array;
%         temp_array(:,simInfo.cw_index) = 1e-10;
%         temp_array(:,simInfo.a_index) = 1e10;

%         valCellCSFour = mat2cell(param_array(:,fourFlags),size(param_array,1),ones(1,size(param_array(:,fourFlags),2)));    
%         valCellCSFour = mat2cell(temp_array,size(temp_array,1),ones(1,size(temp_array,2)));    
        
        % right cycle 4 state        
%         ProductionRateRight = productionRateFunctionFourState(valCellCSFour{:});
%         ProductionRateRight = productionRateFunction(valCellCSFour{:});
%         SharpnessRight = sharpnessFunctionFourState(valCellCSFour{:})*simInfo.crs;
%         SharpnessRight = sharpnessFunction(valCellCSFour{:})*simInfo.crs;        
%         SharpnessRightNorm = SharpnessRight ./ (ProductionRateRight.*(1-ProductionRateRight));
       
        % wrong cycle
%         fourWrongFlags = simInfo.fourStateWrongInputFlags;
%         valCellCSFourWrong = mat2cell(param_array(:,fourWrongFlags),size(param_array,1),ones(1,size(param_array(:,fourWrongFlags),2)));    
%         ProductionRateWrong = productionRateWrongFunctionFourState(valCellCSFourWrong{:});
        
        % inter-network
%         specFactorAlt = log10((ProductionRateRight./ProductionRateWrong)./(simInfo.specFactor/simInfo.cw));
%         deviationFactor = sharpnessWrongFunction(paramCell{:});
%         deviationFactor = 1./deviationFactor;
        
        % intra-network             
        right_flags = simInfo.n_right_bound.*simInfo.activeStateFilter;
        wrong_flags = simInfo.n_wrong_bound.*simInfo.activeStateFilter;
        right_weights = sum(right_flags.*stateProbs,2);
        wrong_weights = sum(wrong_flags.*stateProbs,2);
        specificityFactor = log10((right_weights./wrong_weights)./(simInfo.specFactor./cwVec));                              
        fe_drop = NaN(size(specificityFactor));
        
        % try something
        specF = right_weights./wrong_weights .*cwVec;
        prefactor = specF ./ (specF + cwVec);
        SharpnessRight = Sharpness./prefactor;
%         SharpnessRightAltNorm = SharpnessRightAlt
        
        
    else
        % Free energy drop
        fe_drop = log(prod(param_array(:,simInfo.forwardRateIndices{1}),2)./prod(param_array(:,simInfo.backwardRateIndices{1}),2));
        
        % Cycle rate J                       
%         ProductionRateRight = NaN(size(ProductionRate));
        specificityFactor = NaN(size(ProductionRate));
%         specFactorAlt = NaN(size(ProductionRate));
%         SharpnessRight = specificityFactor; 
%         SharpnessRightNorm = specificityFactor;
%         deviationFactor = specificityFactor;
        cwVec = specificityFactor;
        SharpnessRight = specificityFactor;
%         rate_array_norm = param_array(:,2:end)./ sum(param_array(:,2:end),2);
%         rate_entropy_vec = sum(rate_array_norm.*log(rate_array_norm),2);
%         rate_entropy_vec = rate_entropy_vec - min(rate_entropy_vec);
%         rate_entropy_vec = rate_entropy_vec/nanmax(rate_entropy_vec);
    end
    
    % calculate decision metrics 
    [V, T, R , ~, ~] = calculateDecisionMetrics(param_array,simInfo);                       
    
    % calculate variance that includes poisson noise
    init_rate = simInfo.r_target./ProductionRate/simInfo.TauCycleTime;
    VariancePoisson = Variance.*init_rate.^2 + ProductionRate.*init_rate;
    IRPoisson = Sharpness.^2 ./ VariancePoisson;
    
    % "Affinity" vec
    AffinityVec = log10(TauOff./TauOn);
    
    % calculate entropy rate
    Phi = entropyRateFunction(paramCell{:}).*TauCycle;    
    
    % simple binomial variance
    BinomialVariance = ProductionRate.*(1-ProductionRate);
    
    % generate vector to output
    PrecisionOut = log(BinomialVariance.^2.*TauCycle./Variance);      
    PrecisionRawOut = log(TauCycle./Variance);      
      
    metric_vec = [ProductionRate, Sharpness, SharpnessRaw, fe_drop, PrecisionOut,PrecisionRawOut,...                
                  Phi, TauCycle,specificityFactor,...                                              
                  V.*TauCycle,1./T,R,...                  
                  TauOn,TauOff,AffinityVec,log10(cwVec),...
                  stateEntropyVec,rate_entropy_array,rand_vec,p3Right,p3Wrong,log(TauCycle./VariancePoisson),IRPoisson.*TauCycle,SharpnessRight];
        
    % Apply QC controls
    im_flags =  any(round(real(metric_vec),4)~=round(metric_vec,4) & ~isnan(metric_vec),2);
    val_flags = ProductionRate<=1e-2 | ProductionRate>=(1-1e-2);
    metric_vec(im_flags|val_flags,:) = NaN;
    metric_vec = real(metric_vec);
       
%     if simInfo.saturationFlag==1
%         metric_vec(round(ProductionRateHigh,3)<0.99,:) = NaN;
%     end
else
    metric_vec = [];  
end

metric_names = {'Production Rate','Sharpness','SharpnessRaw','Flux','Precision','PrecisionRaw',...
                'Phi', 'CycleTime','Specificity',...                
                'DecisionRateNorm','DecisionTimeNorm','ProductionBinary',...
                'TauOn','TauOff','AffinityVec','CW','stateEntropy','rateEntropy','rand_vec','p3Right','p3Wrong',...
                'PrecisionPoisson','DecisionRatePoisson','SharpnessRight'};

% specify default bounds
metric_ub_vec = repelem(Inf,length(metric_names));
metric_lb_vec = repelem(-Inf,length(metric_names));
