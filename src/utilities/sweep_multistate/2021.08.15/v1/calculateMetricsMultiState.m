function [metric_vec, metric_names, metric_ub_vec, metric_lb_vec] = ...
                              calculateMetricsMultiState(rate_array, c_val, simInfo)                            
                            
if ~isempty(rate_array)       
  
    % solve for c value where production rate is at half max    
    c1 = c_val + 0.1;
    minError = .32; 
    
    % generate input cell 
    valMat = [c_val*ones(size(rate_array,1),1) rate_array];
    valCell = mat2cell(valMat,size(valMat,1),ones(1,size(valMat,2)));
    
    % calculate production rate         
    ProductionRate = productionRateFunction(valCell{:});
    
    % calculate sharpness     
    Sharpness = sharpnessFunction(valCell{:})*simInfo.c_val;
          
    % calculate variance    
    Variance = intrinsicVarianceFunction(valCell{:});
    
    % get state probabilities
    state_probs = steadyStateVecFunction(valCell{:});
    
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
                
        valMatRight = [c_val*ones(size(rate_array,1),1) rate_array_right];
        valCellRight2 = mat2cell(valMatRight,size(valMatRight,1),ones(1,size(valMatRight,2)));
        
%         SharpnessRightAlt = sharpnessFunction(valCellRight2{:});
%           
%         
        
        % right cycle 4 state
        rightCycleRates = rate_array(:,simInfo.rightIndices);
        valMatRight = [c_val*ones(size(rate_array,1),1) rightCycleRates];
        valCellRight = mat2cell(valMatRight,size(valMatRight,1),ones(1,size(valMatRight,2)));
        ProductionRateRight = productionRateFourStateOR(valCellRight{:});
        SharpnessRight = sharpnessFunctionFourStateOR(valCellRight{:})*simInfo.c_val;%sharpnessRightFunction(valCellRight{:});      
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
    [V,T,R] = calculateDecisionMetrics(rate_array,c1,c_val,minError);
    
    % calculate cycle time
    valCellTrunc = valCell;
%     if simInfo.nStates == 4
%         valCellTrunc = valCellTrunc(:,[5 9]);        
%     end
    TauOn = TauONFunction(valCellTrunc{:});
    TauOff = TauOFFFunction(valCell{:});
    TauCycle = TauOff+TauOn;        
    
    Phi = entropyRateFunction(valCell{:});%.*TauCycle;
%     Phi = CycleFlux.*fe_drop;
    % calculate the ''normalized'' sharpnes as a function of 
%     SharpnessNormed = Sharpness ./ (ProductionRate.*(1-ProductionRate));
    
    valMat1 = [c1*ones(size(rate_array,1),1) rate_array];
    valCell1 = mat2cell(valMat1,size(valMat1,1),ones(1,size(valMat1,2)));   
    
    BinomialVariance = ProductionRate.*(1-ProductionRate);
    % calculate effective non-eq noise reduction
%     V_eq_max = 0.0024;
%     NeqNoiseFactor = V.*tauCycleMean ./ V_eq_max ./ SharpnessNormed.^2;
    
%     SwitchingPrecision = log(tauCycle.*ProductionRate.*(1-ProductionRate)./Variance);
    
    % generate vector to output
    metric_vec = [ProductionRate, Sharpness, fe_drop, log(sqrt(TauCycle./Variance)),...                
                  Phi.*TauCycle, TauCycle,specFactor,ProductionRateRight,SharpnessFactor,PrecisionFactor,...                            
                  SharpnessRight, log(sqrt(TauCycleRight./VarRight)), SharpnessRightNorm, ...
                  log(TauOn./TauOff),V.*TauCycle,1./T,R,...                  
                  log(1./(BinomialVariance)),TauOn,TauOff,specFactorAlt,AffinityVec];
        
    im_flags =  any(round(real(metric_vec),4)~=round(metric_vec,4) & ~isnan(metric_vec),2);
    metric_vec(im_flags,:) = NaN;
    metric_vec = real(metric_vec);
%     param_array(im_flags,:) = NaN;
%     metric_vec(TauCycle<100,:) = NaN;
%     rate_flags1 = abs(log10(rate_array(:,1)./rate_array(:,8))) > 2;
%     metric_vec(rate_flags1,:) = NaN;
%     rate_flags2 = rate_array(:,3)>1e2 | rate_array(:,6)>1e2;
%     metric_vec(rate_flags2,:) = NaN;
else
    metric_vec = [];  
end

metric_names = {'Production Rate','Sharpness','Flux','Precision',...
                'Phi', 'CycleTime','Specificity','ProductionRateRight', 'SharpnessFactor', 'PrecisionFactor',...
                'SharpnessRight', 'PrecisionRight','SharpnessRightNorm','ONOFFRatio',...
                'DecisionRateNorm','DecisionTimeNorm','ProductionBinary',...
                'BinomialNoise','TauOn','TauOff','specFactorAlt','AffinityVec'};

% specify default bounds
metric_ub_vec = repelem(Inf,length(metric_names));
metric_lb_vec = repelem(-Inf,length(metric_names));
