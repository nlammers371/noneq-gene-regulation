function new_rates = applyHMConstraintSixState(simInfo,rate_array,rate_bounds)

% assign other variable values
c = simInfo.c_val;
beta = simInfo.specFactor;
cw = simInfo.wrongFactorConcentration;

% extract subset of unique rate values
rightFlags = ~ismember(1:size(rate_bounds,2),simInfo.rightWrongRates(:,2)');
rightRates = rate_array(:,rightFlags);

inputValMat = [beta*ones(size(rate_array,1),1) c*ones(size(rate_array,1),1)...
                cw*ones(size(rate_array,1),1) rightRates];
% inputValCell = mat2cell(inputValMat,size(inputValMat,1),ones(1,size(inputValMat,2)));

% extract index array
indexArray = simInfo.indArray;
index_cell = {[indexArray(2,1),indexArray(3,4),indexArray(5,4),indexArray(6,1)],...
              [indexArray(1,2),indexArray(4,3),indexArray(4,5),indexArray(1,6)],...
              [indexArray(1,4),indexArray(2,3),indexArray(6,5)],...
              [indexArray(4,1),indexArray(3,2),indexArray(5,6)]};
      
% flags to track which rates still need updating
updateIndices = 1:size(rate_array,1);%true(size(k21));

% initialize array to store update HM rates
new_rates = NaN(size(rate_array));

% randomly select order in which to apply constraints
option_vec = [1:2 4:12];
calc_order = randsample(option_vec,length(option_vec),false);

for i = 1:length(calc_order)
    % get indices
    row = ceil(calc_order(i)/3);
    col = calc_order(i) - 3*(row-1);
    rate_indices = index_cell{row};
    
    % call HM function
    inputValMatIter = inputValMat(updateIndices,:);
    inputValCell = mat2cell(inputValMatIter,size(inputValMatIter,1),ones(1,size(inputValMatIter,2)));
    numStr = [num2str(row) num2str(col)];               
    eval(['aFactor = a' numStr 'SymFun(inputValCell{:});'])
    
    % check validity of solution
    validIndices = find(round(real(aFactor),6) == round(aFactor,6) & real(aFactor)>0);
    
    % now, for subset that qualify, check that they actually lead to HM
    % expression
    rates_temp = rate_array(updateIndices(validIndices),:);
    rates_temp(:,rate_indices) = rates_temp(:,rate_indices).*real(aFactor(validIndices));
    
    tempValMat = [c*ones(size(rates_temp,1),1) rates_temp];
    tempValCell = mat2cell(tempValMat,size(tempValMat,1),ones(1,size(tempValMat,2)));
    
    pdRates = productionRateFunction(tempValCell{:});
    valFlags = round(pdRates,6)== round(1/2,6);
    
    % finally check that rate bounds are obeyed 
    err_flags = max(log10(rates_temp)<rate_bounds(1,:)...
            |log10(rates_temp)>rate_bounds(2,:),[],2);
          
    % update (apologies for the nested indexing)
    new_rates(updateIndices(validIndices(~err_flags&valFlags)),:) = rates_temp(~err_flags&valFlags,:);
    updateIndices = find(isnan(new_rates(:,1)));
end




