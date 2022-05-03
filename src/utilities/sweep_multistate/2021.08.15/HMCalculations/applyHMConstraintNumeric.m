function new_param_array = applyHMConstraintNumeric(param_array, simInfo)

    % set optimization options
    options = optimoptions(@lsqnonlin,'Display','off','MaxIterations',5);
    numerical_precision = sqrt(simInfo.numerical_precision); % let's be alittle more lax for the HM stuff
    % extract param bound info
    paramBounds = simInfo.paramBounds;
   
    % initialize array to store results
    new_param_array = NaN(size(param_array));
    n_trials = size(param_array,1);
    
    % extract elligible rate pairs to sweep
    hmPairs = simInfo.hmRatePairIndices;
    hmPairs = hmPairs(all(simInfo.sweepFlags(hmPairs),2),:);
    
    % randomly select rate pair to adjust for each row
    hmIndices = randsample(1:size(hmPairs,1),n_trials,true);

    for i = 1:n_trials
        % initialize random vector of values
        paramVec = param_array(i,:);        
        if ~any(isnan(paramVec))
            % RAlt = subHelper(valVec,a,10)
            subIndices = unique(hmPairs(hmIndices(i),:));
            add_factor_fun = @(a) subHelper(paramVec,a,subIndices);

            %%define objective function
            ob_fun = @(a) hmCalcHelper(add_factor_fun(a),simInfo.activeStateFilter==1);    
            [a_fit, ~, resid] = lsqnonlin(ob_fun,1,0,1e4,options);

            % update vec
            paramVecNew = paramVec;
            paramVecNew(subIndices) = paramVecNew(subIndices)*a_fit;

            err_flag = max(log10(paramVecNew(simInfo.sweepFlags))<paramBounds(1,simInfo.sweepFlags)...
                        |log10(paramVecNew(simInfo.sweepFlags))>paramBounds(2,simInfo.sweepFlags),[],2);

            if ~err_flag && abs(resid)<=10^-numerical_precision
%                 valCellCS = mat2cell(paramVecNew,size(paramVecNew,1),ones(1,size(paramVecNew,2)));             
% 
%                 % get rate arrays
%                 Q_num_cs = RSymFun(valCellCS{:});
% 
%                 % calculate probabilities
%                 state_probs = calculate_ss_num(Q_num_cs,4);
%                 pdRate = sum(state_probs(simInfo.activeStateFilter==1));
                new_param_array(i,:) = paramVecNew;
            end
        end
    end

    