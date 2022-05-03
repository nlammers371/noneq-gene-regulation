function ssArray = calculateSSProbsNum(simInfo,rate_array,c_val)

    % Initialize
    nStates = simInfo.nStates;
    ssArray = NaN(size(rate_array,1),nStates);
    
    % Solve numerically one by one
%     tic
    for i = 1:size(rate_array,1)
      
        % convert to cell
        valMat = [c_val rate_array(i,:)];
        valCell = mat2cell(valMat,1,ones(1,size(valMat,2)));
        
        % get transition rate matrix
        R_temp = RSymFun(valCell{:});
        
        % calculate eigenvalues/vectors
        [V,D] = eig(R_temp);
        DLog = logical(round(real(D),1)==0);
        ssInd = find(all(DLog));
        ssVec = V(:,ssInd) / sum(V(:,ssInd));
        if all(real(ssVec)==ssVec)
            ssArray(i,:) = ssVec;
        end
    end
%     toc