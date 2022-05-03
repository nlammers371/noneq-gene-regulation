function new_params = applyHMConstraintSymbolic(param_array,simInfo)

hmStruct = simInfo.hmStruct;
paramBounds = simInfo.paramBounds;
sweepFlags = simInfo.sweepFlags;

% get cycle rate indices 
numerical_precision = simInfo.numerical_precision;      

% generate inputy cell
valCell = mat2cell(param_array,size(param_array,1),ones(1,size(param_array,2)));  

%% Make sure we're only adjusting params that are elligible
hm_indices = [];
for i = 1:length(hmStruct)
    if all(simInfo.sweepFlags(hmStruct(i).output_indices))
        hm_indices(end+1) = i;
    end
end

%%
% iterate randomly through different solutions until we find one that works
if simInfo.nStates == 6
    option_vec = [1:2 4:5 7:12];
    calc_order = randsample(option_vec,length(option_vec),false);
elseif simInfo.nStates == 4
    option_vec = [1:8];
    calc_order = randsample(option_vec,length(option_vec),false);
end


% Initialize new rate array
new_params = NaN(size(param_array));

% keep track of indices to update
updateIndices = 1:size(param_array,1);

for i = 1:length(calc_order)
  
    % generate truncated array to update
     inputValMatIter = param_array(updateIndices,:);
     inputValCell = mat2cell(inputValMatIter,size(inputValMatIter,1),ones(1,size(inputValMatIter,2)));
     
    if simInfo.nStates == 4
        % get indices
        row = ceil(calc_order(i)/2);
        col = calc_order(i) - 2*(row-1);        
        % get indices of rates to update
        rate_indices = hmStruct(row).output_indices;
    
        eval(['aFactor = hmFun' num2str(row) num2str(col) '(inputValCell{:});'])
    elseif simInfo.nStates == 6
        % get indices
        row = ceil(calc_order(i)/3);
        col = calc_order(i) - 3*(row-1); 
        
        % get indices of rates to update
        rate_indices = hmStruct(row).output_indices;
        
        % call HM function       
        numStr = [num2str(row) num2str(col)];               
        eval(['aFactor = a' numStr 'SymFun(inputValCell{:});'])
    end
    
    
    % check validity of solution
    validIndices = find(round(real(aFactor),numerical_precision) == round(aFactor,numerical_precision) & real(aFactor)>0);
    
    if ~isempty(validIndices)
        % now, for subset that qualify, check that they actually lead to HM
        % expression
        params_temp = param_array(updateIndices(validIndices),:);       
        params_temp(:,rate_indices) = params_temp(:,rate_indices).*real(aFactor(validIndices));
        
        tempValCell = mat2cell(params_temp,size(params_temp,1),ones(1,size(params_temp,2)));

        pdRates = productionRateFunction(tempValCell{:});
        valFlags = round(pdRates,numerical_precision)== round(1/2,6);

        % finally check that rate bounds are obeyed 
        err_flags = max(log10(params_temp(:,sweepFlags))<paramBounds(1,sweepFlags)...
                |log10(params_temp(:,sweepFlags))>paramBounds(2,sweepFlags),[],2);


        % update (apologies for the nested indexing)
        new_params(updateIndices(validIndices(~err_flags&valFlags)),:) = params_temp(~err_flags&valFlags,:);
        updateIndices = find(isnan(new_params(:,1)));
    end
end

% 
% 
% while ~exit_flag && iter < length(hmStruct)
%   
%     ind = attempt_order(iter);
%     
%     for j = 1:2
%                    
% %         eval(['conditionFlags = conditionsFun' num2str(ind) num2str(j) '(valCell{:});'])    
%         
%         % check which solutions make sense
%         validFlags = real(aValVec) > 0 & round(real(aValVec),numerical_precision)==round(aValVec,numerical_precision) & updateFlags;
%         validIndices = find(validFlags);
%         
%         % apply to rates        
%         update_indices = hmStruct(ind).output_indices;
%         params_temp = param_array(validFlags,:);        
%         params_temp(:,update_indices) = params_temp(:,update_indices).*aValVec(validFlags);
%         
%         % check to see that all new rates are in fact 1/2    
%         valCellTemp = mat2cell(params_temp,size(params_temp,1),ones(1,size(params_temp,2))); 
%         pdRateTemp = productionRateFunction(valCellTemp{:});    
%         hasCorrectValue = round(pdRateTemp,numerical_precision)==round(0.5,numerical_precision);
%         
%         % check to see if rates violoate magnitude bounds
%         err_flags = max(log10(params_temp)<paramBounds(1,:)...
%                                 |log10(params_temp)>paramBounds(2,:),[],2);
%                     
%         
%         % update factor and index sets 
%         a_factor_vec(validIndices(hasCorrectValue&~err_flags)) = aValVec(validIndices(hasCorrectValue&~err_flags));
% 
%         new_params(validIndices(hasCorrectValue&~err_flags),:) = params_temp(hasCorrectValue&~err_flags,:);
%         
%         % update the update flags
%         updateFlags = isnan(new_params(:,1));
%         
%         exit_flag = all(~updateFlags);
%                 
%     end
%     
%     iter = iter + 1;            
%     
% end
% 
% 
% % instances where we did not find correct solution
% error_flags_val = updateFlags;%round(ProductionRate,5)~=0.500;
% 
% % instances where solution required out-of-bound values
% err_flags_bound = max(log10(new_params)<paramBounds(1,:)...
%                       |log10(new_params)>paramBounds(2,:),[],2);

% replace problematic rates with NaNs
% new_params(error_flags_val|err_flags_bound,:) = NaN;
