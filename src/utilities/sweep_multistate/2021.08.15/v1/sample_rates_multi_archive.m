function new_rates = sample_rates_multi_archive(lb_array,ub_array,orig_rate_array,prop_sigma,simInfo)

  % extract key parameters
  rate_bounds = simInfo.rate_bounds;
  nStates = simInfo.nStates; 
  equilibrium_flag = simInfo.equilibrium_flag;
  
  % generate proposed rates
  if isempty(orig_rate_array)
      new_rates = reshape(10.^(prop_sigma*trandn(lb_array,ub_array)),[],3*nStates-4); % function to sample from truncated normal dist
  else
      new_rates = 10.^(orig_rate_array+reshape(prop_sigma*trandn(lb_array,ub_array),[],3*nStates-4));
  end   
        
  % apply symmetry constraints
  if simInfo.nStates > 4
      new_rates(:,simInfo.rightWrongRates(:,2)) = new_rates(:,simInfo.rightWrongRates(:,1));
  end
  
  % check to see if there are any special constraints to apply
  if any(simInfo.unbindingFlagsWrong)
      ubWrongMatchArray = simInfo.ubWrongMatchArray;
      new_rates(:,ubWrongMatchArray(:,2)) = new_rates(:,ubWrongMatchArray(:,1))*simInfo.specFactor;      
      new_rates(:,simInfo.bindingFlagsWrong) = new_rates(:,simInfo.bindingFlagsWrong)*simInfo.wrongFactorConcentration;
  end
    
  if equilibrium_flag
    
      % generate input cell 
      c_val = 1;
      valMat = [c_val*ones(size(new_rates,1),1) new_rates];
      valCell = mat2cell(valMat,size(valMat,1),ones(1,size(valMat,2)));
      
      % calculate state probabilities
      state_probs = steadyStateVecFunction(valCell{:});          
      
      % extract useful indexing arrays       
      connectedStatePairs = simInfo.connectedStatePairs;
      fromRateArray = simInfo.fromRateArray;
      
      % use rejection sampling to find reverse rates
      new_rates_temp = new_rates;
      err_flags = true(size(state_probs,1),1);
      iter = 0; 
      validFraction = 0;
      
      while iter < 10 && validFraction < 0.95             

          % attempt to resample reverse rates         
          state_prob_ratios = state_probs(err_flags,connectedStatePairs(:,1))./state_probs(err_flags,connectedStatePairs(:,2));
          new_rates_temp(err_flags,fromRateArray(:,2)) = state_prob_ratios.*new_rates_temp(err_flags,fromRateArray(:,1));
  
          % find cases where this renormalization pushes rates outside of bounds
          err_flags = max(log10(new_rates_temp(:,fromRateArray(:,2)))<rate_bounds(1,fromRateArray(:,2))...
            |log10(new_rates_temp(:,fromRateArray(:,2)))>rate_bounds(2,fromRateArray(:,2)),[],2);
          
          % resample state probabilities for cases where rates fall outside
          % of bounds
          state_probs(err_flags,:) = rand(sum(err_flags),nStates);
          state_probs(err_flags,:) = state_probs(err_flags,:) ./ sum(state_probs(err_flags,:),2);
          
          % calculate metrics
          iter = iter + 1;
          validFraction = mean(~err_flags);
      end

%     % reset
      new_rates = new_rates_temp;
      if simInfo.nStates > 4
          new_rates(:,simInfo.rightWrongRates(:,2)) = new_rates(:,simInfo.rightWrongRates(:,1));
      end
      if any(simInfo.unbindingFlagsWrong)
          ubWrongMatchArray = simInfo.ubWrongMatchArray;
          new_rates(:,ubWrongMatchArray(:,2)) = new_rates(:,ubWrongMatchArray(:,1))*simInfo.specFactor;
          new_rates(:,simInfo.bindingFlagsWrong) = new_rates(:,simInfo.bindingFlagsWrong)*simInfo.wrongFactorConcentration;
      end
      
      new_rates(err_flags,:) = NaN;
    
  end
  

  
  
  
