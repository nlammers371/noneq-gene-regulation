function new_rates = applyHMConstraint(simInfo,rate_array,rate_bounds)

hmStruct = simInfo.hmStruct;

% get cycle rate indices 
c_val = simInfo.c_val;
      
% generate inputy cell
valMat = [c_val*ones(size(rate_array,1),1) rate_array];
valCell = mat2cell(valMat,size(valMat,1),ones(1,size(valMat,2)));    

% iterate randomly through different solutions until we find one that works
attempt_order = randperm(length(hmStruct),length(hmStruct));
exit_flag = false;
updateFlags = true(size(rate_array,1),1);
iter = 1;

% Initialize new rate array
new_rates = rate_array;
% Initialize array to stor multiplier values
a_factor_vec = NaN(size(rate_array,1),1);
rate_index_array = NaN(size(rate_array,1),2);

while ~exit_flag && iter < length(hmStruct)
  
    ind = attempt_order(iter);
    
    % check to see if conditions for solution are met
%     cdFun = hmStruct(ind).conditions;
%     eval(['conditionFlags = conditionsFun' num2str(ind) '1(valCell{:});'   ])
        
    % extract function and solve
%     hmFun = hmStruct(ind).hmSolutions;
%     eval(['aValVec = hmFun' num2str(ind) '1(valCell{:});'   ])
%     aValVec = reshape(hmFun(valCell{:}),[],2);        
    
    for j = 1:2
                   
        eval(['conditionFlags = conditionsFun' num2str(ind) num2str(j) '(valCell{:});'])       
        eval(['aValVec = hmFun' num2str(ind) num2str(j) '(valCell{:});'])
        
        % check which solutions make sense
        validFlags = conditionFlags & real(aValVec) > 0 & real(aValVec)==aValVec;
    
        % apply to rates        
        update_indices = hmStruct(ind).output_indices;
        new_rates(updateFlags&validFlags,update_indices) = ...
              new_rates(updateFlags&validFlags,update_indices).*aValVec(updateFlags&validFlags);
        
            
        % check to see if rates violoate magnitude bounds
        err_flags = max(log10(new_rates(:,update_indices))<rate_bounds(1,update_indices)...
            |log10(new_rates(:,update_indices))>rate_bounds(2,update_indices),[],2);
        
        % update factor and index sets 
        a_factor_vec(updateFlags&validFlags&~err_flags) = aValVec(updateFlags&validFlags&~err_flags);
        rate_index_array(updateFlags&validFlags&~err_flags,:) = repmat(update_indices,sum(updateFlags&validFlags&~err_flags),1);
        
        % update the update flags
        updateFlags = (updateFlags & ~validFlags) | err_flags;
        
        exit_flag = all(~updateFlags);
                
    end
    
    iter = iter + 1;            
    
end

% check for errors
% valMatNew = [c_val*ones(size(rate_array,1),1) new_rates];
% valCellNew = mat2cell(valMatNew,size(valMatNew,1),ones(1,size(valMatNew,2)));
% ProductionRate = productionRateFunction(valCellNew{:});

% instances where we did not find correct solution
error_flags_val = updateFlags;%round(ProductionRate,5)~=0.500;

% instances where solution required out-of-bound values
err_flags_bound = max(log10(new_rates)<rate_bounds(1,:)...
                      |log10(new_rates)>rate_bounds(2,:),[],2);

% replace problematic rates with NaNs
new_rates(error_flags_val|err_flags_bound,:) = NaN;
% a_factor_vec(error_flags_val|err_flags_bound) = NaN;
% rate_index_array(error_flags_val|err_flags_bound,:) = 1; % to prevent indexing error