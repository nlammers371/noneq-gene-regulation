function rates_out = constrainCycleTime(rate_array, c_range, ...
                                                            rate_bounds, T)

for k = 1:4
    eval(['k' num2str(k) ' = rate_array(:,k);'])
end
for r = 1:4
    eval(['r' num2str(r) ' = rate_array(:,r+4);'])
end

% get indices of edges that are elligible for adjustment
rates_out = rate_array;

replacementFactors = calculateCycleRates(rate_array, c_range, T);
replacementRates = rates_out.*replacementFactors;

% look for proposed values that meet criteria
filterArray = isreal(replacementRates) & replacementRates >= 10.^rate_bounds(1,:) & replacementRates <= 10.^rate_bounds(2,:);  
filterArray = filterArray(:,1:4) & filterArray(:,5:8);

nOptions = sum(filterArray,2);
choiceMat = cumsum(filterArray,2) ./ nOptions;
[~,decisionVec] = max(choiceMat > rand(size(nOptions)),[],2);
linIndices = sub2ind(size(rate_array),[1:length(decisionVec) 1:length(decisionVec)],[decisionVec' decisionVec'+4]);

rates_out(linIndices) = replacementRates(linIndices);
rates_out(nOptions==0,:) = NaN;
% iterate through each set of rates and randomly select edge to replace
% for i = 1:size(rate_array,1)
%     edgeOptions = find(filterArray(i,:));
%     if isempty(edgeOptions)
%         rates_out(i, :) = NaN;
%     elseif length(edgeOptions)==1
%         rates_out(i, [edgeOptions edgeOptions+4]) = replacementRates(i,[edgeOptions edgeOptions+4]);
%     else
%         ri = randsample(edgeOptions,1);
%         rates_out(i, [ri ri+4]) = replacementRates(i,[ri ri+4]);
%     end
% end
