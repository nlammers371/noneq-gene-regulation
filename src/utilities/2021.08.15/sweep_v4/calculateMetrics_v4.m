function [metric_vec, metric_names,metric_ub_vec,metric_lb_vec] = calculateMetrics_v4(rate_array, c_val)

if ~isempty(rate_array)            
    % solve for c value where production rate is at half max    
    c1 = 1.1*c_val;
    minError = .1; 
    
    % generate input cell 
    valMat = [c_val*ones(size(rate_array,1),1) rate_array];
    valCell = mat2cell(valMat,size(valMat,1),ones(1,size(valMat,2)));
    
    % calculate production rate         
    ProductionRate = productionRateFunction4State(valCell{:});
    
    % calculate sharpness     
    Sharpness = sharpnessFunction4State(valCell{:});
        
    if size(rate_array,2)==8
        % calculate thermodynamic force
        fe_drop = log(prod(rate_array(:,1:4),2)./prod(rate_array(:,5:8),2));
        % calculate unidirectional cycle fluxes    
        Jp = forwardFluxFunction4State(valCell{:});
        Jm = backwardFluxFunction4State(valCell{:});
    else
        fe_drop = NaN;
        Jp = NaN;
        Jm = NaN;
    end
        
    % calculate variance    
    VarPoint = intrinsicVarianceFunction4State(valCell{:});
    
    % calculate decision metrics 
    [V,T,R] = calculateDecisionMetrics(rate_array,c1,c_val,minError);
    
    % calculate cycle time
%     [TauOn,TauOff,~,~,dKondC,dKoffdC] = fourStateCycleTime_v4(rate_array,c_val);            
    
    TauOn = TauONFunction4State(valCell{:});
    TauOff = TauOFFFunction4State(valCell{:});
    tauCycle = TauOff+TauOn;
    
    % Calculate the average fraction of time the system spends in each
    % state    
    steadyStateVec = reshape(steadyStateVecFunction4State(valCell{:}),[],4);
    
    % calculate the ''normalized'' sharpnes as a function of 
    SharpnessNormed = Sharpness ./ (ProductionRate.*(1-ProductionRate));
    
    valMat1 = [c1*ones(size(rate_array,1),1) rate_array];
    valCell1 = mat2cell(valMat1,size(valMat1,1),ones(1,size(valMat1,2)));
    
    TauOn1 = TauONFunction4State(valCell1{:});
    TauOff1 = TauOFFFunction4State(valCell1{:});     
    tauCycle1 = TauOff1+TauOn1;
    tauCycleMean = 0.5*tauCycle + 0.5*tauCycle1;
    
    % calculate effective non-eq noise reduction
    V_eq_max = 0.0024;
    NeqNoiseFactor = V.*tauCycleMean ./ V_eq_max ./ SharpnessNormed.^2;
    
    % generate vector to output
    metric_vec = [ProductionRate, Sharpness, fe_drop, log(sqrt(tauCycle./VarPoint))...                
                  Jp, Jm, (Jp-Jm).*fe_drop, Jp+Jm, tauCycle,...                  
                  steadyStateVec(:,1)+steadyStateVec(:,3),...
                  SharpnessNormed, log(TauOn./TauOff),V.*tauCycleMean,tauCycleMean./T,R,V,...
                  log(tauCycle.*ProductionRate.*(1-ProductionRate)./VarPoint),...
                  log(1./(ProductionRate.*(1-ProductionRate))),TauOn,TauOff,NeqNoiseFactor];
        
else
    metric_vec = [];  
end

metric_names = {'Production Rate','Sharpness','Flux','Precision',...
                'Positive Flux','Negative Flux', 'Phi', '|J|', 'CycleTime',...
                '2StateFrac','SharpnessNormed','Affinity',...
                'DecisionRateNorm','DecisionTimeNorm','ProductionBinary',...
                'DecisionRate','SwitchingNoise','BinomialNoise','TauOn','TauOff','NeqNoiseFactor'};

% specify default bounds
metric_ub_vec = repelem(Inf,length(metric_names));
metric_lb_vec = repelem(-Inf,length(metric_names));
