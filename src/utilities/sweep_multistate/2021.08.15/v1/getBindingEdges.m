function simInfo = getBindingEdges(simInfo)   

    nStates = simInfo.nStates;
    % this indexing function assumes the presence of a mirrored "wrong"
    % half of the graph, so for the simple 4 state binding case we include
    % a shadow "wrong" side and just zero out the connection matrices 
    % accordingly 
    fourFlag = 0;    
    if nStates == 4
        fourFlag = 1;
        nStates = 6;
    end
    
    % NL: hard-coding right forward and back rates for now 
    if fourFlag
        simInfo.forwardRatesRight = [3 5 8 2];
        simInfo.backwardRatesRight = [7 6 4 1];
    else
        simInfo.forwardRatesRight = [4 6 9 2];
        simInfo.backwardRatesRight = [8 7 5 1];
    end
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
    permittedConnections = permittedConnections(si,si);
    if fourFlag
        permittedConnections(5:end,:) = 0;
        permittedConnections(:,5:end) = 0;
    end
    rateIndices= find(permittedConnections)';
    [~,siFunction] = sort(toRef(permittedConnections));
    rateIndices = rateIndices(siFunction);

    bindingConnections = bindingConnections(si,si);
    unbindingConnections = unbindingConnections(si,si);
    if fourFlag
        unbindingConnections(5:end,:) = 0;
        unbindingConnections(:,5:end) = 0;
        bindingConnections(5:end,:) = 0;
        bindingConnections(:,5:end) = 0;
    end
    unbindingConnectionsWrong = unbindingConnections & ~fromRefHalved;
    bindingConnectionsWrong = bindingConnections & ~toRefHalved;

    simInfo.bindingFlags = bindingConnections(rateIndices);
    simInfo.unbindingFlags = unbindingConnections(rateIndices);
    simInfo.bindingFlagsWrong = bindingConnectionsWrong(rateIndices);
    simInfo.unbindingFlagsWrong = unbindingConnectionsWrong(rateIndices);
    ubWrongIndices = find(simInfo.unbindingFlagsWrong);
    
    % get matching states
    vec = 1:nStates;
    distVec = min([abs(vec'-vec(1)) abs(vec'-vec(end)-1)],[],2)';
    fromVec = fromRef(rateIndices);
    ubWrongDistances = distVec(fromVec(simInfo.unbindingFlagsWrong));
    ubRightMatches = find(ismember(distVec(fromVec),ubWrongDistances)&~simInfo.unbindingFlagsWrong&simInfo.unbindingFlags);
    
    simInfo.ubWrongMatchArray = [ubRightMatches' ubWrongIndices'];
    
    % get connected node info
    A = triu(permittedConnections);
    [toStates, fromStates] = find(A);
    simInfo.connectedStatePairs = unique([toStates, fromStates],'rows');

    % create and index array
    simInfo.indArray = zeros(nStates);
    simInfo.indArray(rateIndices) = 1:length(rateIndices);

    simInfo.fromRateArray = zeros(size(simInfo.connectedStatePairs));
    for c = 1:size(simInfo.connectedStatePairs,1)
        simInfo.fromRateArray(c,1) = simInfo.indArray(simInfo.connectedStatePairs(c,2),simInfo.connectedStatePairs(c,1));
        simInfo.fromRateArray(c,2) = simInfo.indArray(simInfo.connectedStatePairs(c,1),simInfo.connectedStatePairs(c,2));
    end
    
    %% calculate right-wrong matches to enforce symmetries
    if ~fourFlag
        rightVals = distVec(1:nStates/2);
        wrongVals = distVec(nStates/2+1:end);
        stateMatchArray = [];
        for u = 2:length(rightVals)
            stateMatchArray = vertcat(stateMatchArray,[u,nStates/2 + find(wrongVals==rightVals(u))]);
        end

        % cross-reference with connected state array to get right/wrong edge
        % matches
        cpMat = simInfo.connectedStatePairs;
        fMat = simInfo.fromRateArray;
        wrongStates = stateMatchArray(:,2);
        rightStates = stateMatchArray(:,1);
        simInfo.rightWrongRates = NaN(floor(nStates/4)*6,2);
        iter = 1;
        for c = 1:size(cpMat,1)
            row = cpMat(c,:);
            repIndices = find(ismember(row,wrongStates));
            if ~isempty(repIndices)
                rowNew = row;
                for r = repIndices
                   rowNew(r) = rightStates(row(r)==wrongStates);
                end
                [rowNew, si] = sort(rowNew);
                row_ref = all(ismember(cpMat,rowNew),2);
                simInfo.rightWrongRates(iter:iter+1,1) = fMat(row_ref,:);
                simInfo.rightWrongRates(iter:iter+1,2) = fMat(c,si);
                iter = iter + 2;
            end
        end     
        
        % generate flags that map 6 state rates to 4 state networks for
        % "right" and "wrong" cycles
        rightIndices = simInfo.indArray(1:4,1:4)';
        simInfo.rightIndices = rightIndices(rightIndices~=0);
        
        wrongIndices = simInfo.indArray([1 6 5 4],[1 6 5 4])';
        simInfo.wrongIndices = wrongIndices(wrongIndices~=0);
    else       
        simInfo.rightWrongRates = NaN(0,2);
    end
    
    