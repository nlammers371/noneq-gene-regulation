function [metric_vec, metric_names, metric_ub_vec, metric_lb_vec] = ...
                              calculateMetricsMultiStateNumSS(rate_array, c_val, simInfo)                            
                            
if ~isempty(rate_array)       
  
    % initialize parameter arrays
    state_probs_c0 = NaN(size(rate_array,1),simInfo.nStates);  
    state_probs_cs = NaN(size(rate_array,1),simInfo.nStates);  
    state_probs_c1 = NaN(size(rate_array,1),simInfo.nStates);  
    Variance = NaN(size(rate_array,1),1);
    TauCycle = NaN(size(rate_array,1),1);
    TauOn = NaN(size(rate_array,1),1);
    TauOff = NaN(size(rate_array,1),1);
    Phi = NaN(size(rate_array,1),1);
    
    % generate input arrays
    val = rand(1,11);
    valCell = mat2cell(val,size(val,1),ones(1,size(val,2)));
    for i = 1:size(rate_array,1)
        
        activeStateFilter = ismember(1:nStates,simInfo.activeStates);

        Q_num = RSymFun(valCell{:});
        state_probs(i,:) = calculate_ss_num(Q_num);        
    
        Z_num = calculate_Z_matrix(Q_num,ss_num);
        Variance(i) = calculate_var_num(Q_num,ss_num,activeStateFilter,Z_num);
        Phi(i) = calculate_entropy_rate_num(Q_num,ss_num);
        
        [TauOn(i),TauOff(i),TauCycle(i)] = calculate_tau_num(Q_num,ss_num,activeStateFilter);
      
    end    
    
    % calculate sharpness and production rate
    ProductionRate = sum(state_probs(:,activeStateFilter),2);
    state_probs_alt = calculateSSProbsNum(simInfo,rate_array,c_val*1.05);
    ProductionRateAlt = sum(state_probs_alt(:,activeStates),2);
    Sharpness = (ProductionRateAlt-ProductionRate)./(c_val*1.05-c_val);

    % generate input array  
    valMat = [c_val*ones(size(rate_array,1),1) rate_array state_probs];
    valCell = mat2cell(valMat,size(valMat,1),ones(1,size(valMat,2)));
    
    % calculate variance    
    Variance = intrinsicVarianceFunction(valCell{:});

    % generate pseudo four state network to calculate specificity
    if simInfo.nStates > 4
        % hard-code cycle rate indices for now
        % note that the FE drop will be identical for right and wrong cycles
        forward_right_rates = rate_array(:,simInfo.forwardRatesRight);
        backward_right_rates = rate_array(:,simInfo.backwardRatesRight);
        fe_drop = log(prod(forward_right_rates,2)./prod(backward_right_rates,2));
               
        rate_array_right = rate_array;
        rate_array_right(:,simInfo.rightWrongRates(:,2)) = 1e50;
        rate_array_right(:,simInfo.bindingFlagsWrong) = 1e-50;                     

        % right cycle 4 state
        rightCycleRates = rate_array(:,simInfo.rightIndices);
        valMatRight = [c_val*ones(size(rate_array,1),1) rightCycleRates];
        valCellRight = mat2cell(valMatRight,size(valMatRight,1),ones(1,size(valMatRight,2)));
        ProductionRateRight = productionRateFourStateOR(valCellRight{:});
        SharpnessRight = sharpnessFunctionFourStateOR(valCellRight{:});%sharpnessRightFunction(valCellRight{:});      
        SharpnessRightNorm = SharpnessRight ./ (ProductionRateRight.*(1-ProductionRateRight));
        VarRight = intrinsicVarianceFourStateOR(valCellRight{:});
        bin_factor = (ProductionRateRight.*(1-ProductionRateRight));
        VarRight = VarRight./ bin_factor.^2;
        
        % right cycle 4 state
        wrongCycleRates = rate_array(:,simInfo.wrongIndices);
        valMatWrong = [c_val*ones(size(rate_array,1),1) wrongCycleRates];
        valCellWrong = mat2cell(valMatWrong,size(valMatWrong,1),ones(1,size(valMatWrong,2)));
        ProductionRateWrong = productionRateFunctionFourStateAnd(valCellWrong{:});
       
        % calculate effective affinity for 4 state right cycle 
        onDwellTime = TauONFunctionBinding(valCellRight{:});
        offDwellTime = TauOFFFunctionBinding(valCellRight{:});
        TauCycleRight = offDwellTime + onDwellTime;
        AffinityVec = log10(offDwellTime./onDwellTime);
%         AffinityVec = log10(rightCycleRates(:,1)./rightCycleRates(:,3));
        
        specFactor = log10((state_probs(:,3)./state_probs(:,5))./(simInfo.specFactor/simInfo.wrongFactorConcentration));
        SharpnessFactor = SharpnessRight;%state_probs(:,3)./(state_probs(:,3)+state_probs(:,5)).*SharpnessRight;
        PrecisionFactor = VarRight./Variance;
%         CycleFlux = NaN(size(PrecisionFactor));
        specFactorAlt = log10((ProductionRateRight./ProductionRateWrong)./(simInfo.specFactor/simInfo.wrongFactorConcentration));
        
    else
        forward_right_rates = rate_array(:,simInfo.forwardRatesRight);
        backward_right_rates = rate_array(:,simInfo.backwardRatesRight);
        fe_drop = log(prod(forward_right_rates,2)./prod(backward_right_rates,2));
      
        ProductionRateRight = NaN(size(ProductionRate));
        specFactor = Inf(size(ProductionRate));
        specFactorAlt = Inf(size(ProductionRate));
        SharpnessRight = specFactor;
        SharpnessFactor = specFactor;
        PrecisionFactor = specFactor;
        VarRight = specFactor;
        CycleFlux = state_probs(:,1).*rate_array(:,3) - state_probs(:,2).*rate_array(:,1);
        
        onDwellTime = specFactor;%TauONFunctionBinding(valCell{:});
        offDwellTime = specFactor;%TauOFFFunctionBinding(valCell{:});
        AffinityVec = log10(offDwellTime./onDwellTime);
        TauCycleRight = specFactor;
        SharpnessRightNorm = specFactor;
    end
    % calculate entropy rate
%     Phi = entropyRateFunction(valCell{:});
    
    % calculate decision metrics 
%     [V,T,R] = calculateDecisionMetrics(rate_array,c1,c_val,minError);
    V = NaN(size(Variance));
    R = V;
    T = V;


    tauCycleMean = TauCycle;
    
    BinomialVariance = ProductionRate.*(1-ProductionRate);
    % calculate effective non-eq noise reduction

    
    % generate vector to output
    metric_vec = [ProductionRate, Sharpness, fe_drop, log(sqrt(TauCycle./Variance)),...                
                  Phi, TauCycle,specFactor,ProductionRateRight,SharpnessFactor,PrecisionFactor,...                            
                  SharpnessRight, log(sqrt(TauCycleRight./VarRight)) , SharpnessRightNorm, ...
                  log(TauOn./TauOff),V.*tauCycleMean,tauCycleMean./T,R,V,...
                  ...
                  log(1./(BinomialVariance)),TauOn,TauOff,specFactorAlt,AffinityVec];
        
else
    metric_vec = [];  
end

metric_names = {'Production Rate','Sharpness','Flux','Precision',...
                'Phi', 'CycleTime','Specificity','ProductionRateRight', 'SharpnessFactor', 'PrecisionFactor',...
                'SharpnessRight', 'PrecisionRight','SharpnessRightNorm','ONOFFRatio',...
                'DecisionRateNorm','DecisionTimeNorm','ProductionBinary',...
                'DecisionRate','BinomialNoise','TauOn','TauOff','specFactorAlt','AffinityVec'};

% specify default bounds
metric_ub_vec = repelem(Inf,length(metric_names));
metric_lb_vec = repelem(-Inf,length(metric_names));
