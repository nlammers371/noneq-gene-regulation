function new_params = applyHMConstraintSymbolic(param_array,simInfo)

hmStruct = simInfo.hmStruct;
paramBounds = simInfo.paramBounds;

% get cycle rate indices 
numerical_precision = simInfo.numerical_precision;      

% generate inputy cell
valCell = mat2cell(param_array,size(param_array,1),ones(1,size(param_array,2)));  

% iterate randomly through different solutions until we find one that works
attempt_order = randperm(length(hmStruct),length(hmStruct));
exit_flag = false;
updateFlags = true(size(param_array,1),1);
iter = 1;

% Initialize new rate array
new_params = NaN(size(param_array));

% Initialize array to store multiplier values
a_factor_vec = NaN(size(param_array,1),1);

while ~exit_flag && iter < length(hmStruct)
  
    ind = attempt_order(iter);
    
    for j = 1:2
                   
%         eval(['conditionFlags = conditionsFun' num2str(ind) num2str(j) '(valCell{:});'])       
        eval(['aValVec = hmFun' num2str(ind) num2str(j) '(valCell{:});'])
        
        % check which solutions make sense
        validFlags = real(aValVec) > 0 & round(real(aValVec),numerical_precision)==round(aValVec,numerical_precision) & updateFlags;
        validIndices = find(validFlags);
        
        % apply to rates        
        update_indices = hmStruct(ind).output_indices;
        params_temp = param_array(validFlags,:);        
        params_temp(:,update_indices) = params_temp(:,update_indices).*aValVec(validFlags);
        
        % check to see that all new rates are in fact 1/2    
        valCellTemp = mat2cell(params_temp,size(params_temp,1),ones(1,size(params_temp,2))); 
        pdRateTemp = productionRateFunction(valCellTemp{:});    
        hasCorrectValue = round(pdRateTemp,numerical_precision)==round(0.5,numerical_precision);
        
        % check to see if rates violoate magnitude bounds
        err_flags = max(log10(params_temp)<paramBounds(1,:)...
                                |log10(params_temp)>paramBounds(2,:),[],2);
                    
        
        % update factor and index sets 
        a_factor_vec(validIndices(hasCorrectValue&~err_flags)) = aValVec(validIndices(hasCorrectValue&~err_flags));

        new_params(validIndices(hasCorrectValue&~err_flags),:) = params_temp(hasCorrectValue&~err_flags,:);
        
        % update the update flags
        updateFlags = isnan(new_params(:,1));
        
        exit_flag = all(~updateFlags);
                
    end
    
    iter = iter + 1;            
    
end


% instances where we did not find correct solution
error_flags_val = updateFlags;%round(ProductionRate,5)~=0.500;

% instances where solution required out-of-bound values
err_flags_bound = max(log10(new_params)<paramBounds(1,:)...
                      |log10(new_params)>paramBounds(2,:),[],2);

% replace problematic rates with NaNs
new_params(error_flags_val|err_flags_bound,:) = NaN;
