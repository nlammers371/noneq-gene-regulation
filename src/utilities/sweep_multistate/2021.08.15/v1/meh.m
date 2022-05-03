
nStates = 6;

% zero out forbidden transitions
[fromRef,toRef] = meshgrid(1:nStates,1:nStates);
diffRef1 = fromRef-toRef;
diffRef2 = toRef-fromRef;

% delineate top and bottom halves of array
toRefHalved = toRef<=nStates/2;
fromRefHalved = fromRef<=nStates/2;

% quarters
fromRefQuarters = fromRef<nStates/4 & fromRefHalved;
fromRefQuarters = fromRefQuarters | fromRef<3*nStates/4&~fromRefHalved;
toRefQuarters = toRef<nStates/4 & toRefHalved;
toRefQuarters = toRefQuarters | toRef<3*nStates/4&~toRefHalved;

% generate logical matrices
permittedConnections = (abs(diffRef1)==1 & toRefHalved==fromRefHalved) | abs(diffRef1)==nStates/2;
unbindingConnections = (diffRef1==-1 & fromRefQuarters) | (diffRef2==-1 & ~toRefQuarters) & permittedConnections;
bindingConnections = ~unbindingConnections & permittedConnections & ~(abs(diffRef1)==nStates/2);

% permute these connections to follow a more intuitive labeling scheme
indexCurr = 1:nStates;
indexAdjusted = circshift(indexCurr,floor(nStates/4));
indexAdjusted = [indexAdjusted(1:nStates/2) fliplr(indexAdjusted(nStates/2+1:end))];
[~, si] = sort(indexAdjusted);

% convert these to vectors, noting that the functions will sort variables
% alphabetically. This means that rate variables are sorted by row, then
% column...
rateIndices = find(permittedConnections(si,si))';
[~,siFunction] = sort(toRef(permittedConnections(si,si)));
rateIndices = rateIndices(siFunction);

bindingConnections = bindingConnections(si,si);
unbindingConnections = unbindingConnections(si,si);
unbindingConnectionsWrong = unbindingConnections & ~fromRefHalved;

bindingFlags = bindingConnections(rateIndices);
unbindingFlags = unbindingConnections(rateIndices);
unbindingFlagsWrong = unbindingConnectionsWrong(rateIndices);
ubWrongIndices = find(unbindingFlagsWrong);

RSymFull = sym('k%d%d', [nStates nStates],'positive');
% %%
% test = RSymFull(rateIndices)';
% test(bindingFlags)
%%
vec = 1:nStates;
distVec = min([abs(vec'-vec(1)) abs(vec'-vec(end)-1)],[],2)';
fromVec = fromRef(rateIndices);
ubWrongDistances = distVec(fromVec(unbindingFlagsWrong));
ubRightMatches = find(ismember(distVec(fromVec),ubWrongDistances)&~unbindingFlagsWrong&unbindingFlags);


%% 

% get connected node info
A = triu(permittedConnections(si,si));
[toStates,fromStates] = find(A);
simInfo.connectedStatePairs = unique([toStates, fromStates],'rows');

% create and index array
indArray = zeros(nStates);
simInfo.indArray(rateIndices) = 1:length(rateIndices);

simInfo.fromRateArray = zeros(size(simInfo.connectedStatePairs));
for c = 1:size(simInfo.connectedStatePairs,1)
    simInfo.fromRateArray(c,1) = simInfo.indArray(simInfo.connectedStatePairs(c,2),simInfo.connectedStatePairs(c,1));
    simInfo.fromRateArray(c,2) = simInfo.indArray(simInfo.connectedStatePairs(c,1),simInfo.connectedStatePairs(c,2));
end





