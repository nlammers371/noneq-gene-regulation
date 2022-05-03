function new_rates = applyHMConstraintFourStateNum(simInfo,rate_array,rate_bounds)

    rate_pairs = [3 6; 1 8 ; 4 2; 5 7];
    syms a positive
     
%     rates_temp = rate_array;
%     rates_temp(:,rate_pairs(1,:)) = rates_temp(:,rate_pairs(1,:))*a;
    
   
    
    nStates = simInfo.nStates;
    ssArray = NaN(size(rate_array,1),nStates);
    c_val = simInfo.c_val;
    % Solve numerically one by one
%     tic
    for i = 1:size(rate_array,1)
      
        % convert to cell
        valMat = [c_val rate_array(i,:)];
        valCell = mat2cell(valMat,1,ones(1,size(valMat,2)));
        
        % get transition rate matrix
        R_temp = RSymFun(valCell{:});
        R_sym = sym(R_temp);
        R_sym(2,1) = R_sym(2,1)*a;
        R_sym(3,4) = R_sym(3,4)*a;
        
        %
        % calculate eigenvalues/vectors
        [V,D] = eig(R_sym);
        DLog = logical(round(real(D),1)==0);
        ssInd = find(all(DLog));
        ssVec = V(:,ssInd) / sum(V(:,ssInd));
        
        if all(real(ssVec)==ssVec)
            ssArray(i,:) = ssVec;
        end
    end