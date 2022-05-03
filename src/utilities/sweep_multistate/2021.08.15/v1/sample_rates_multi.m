function new_rates = sample_rates_multi(lb_array,ub_array,orig_rate_array,prop_sigma,simInfo)

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
        
%   if simInfo.sameOnRateFlag
% %       new_rates(:,3) = 1;
% %       new_rates(:,6) = new_rates(:,3);
%       new_rates(:,8) = new_rates(:,1);
%   end
  % apply symmetry constraints
  if simInfo.nStates > 4
      new_rates(:,simInfo.rightWrongRates(:,2)) = new_rates(:,simInfo.rightWrongRates(:,1));
  end
%   new_rates(:,simInfo.bindingFlags) = 1;
  % check to see if there are any special constraints to apply
%   new_rates(:,simInfo.bindingFlags) = new_rates(:,simInfo.bindingFlags).*simInfo.c_val; 
  if any(simInfo.unbindingFlagsWrong)
      ubWrongMatchArray = simInfo.ubWrongMatchArray;
      new_rates(:,ubWrongMatchArray(:,2)) = new_rates(:,ubWrongMatchArray(:,1))*simInfo.specFactor;      
      new_rates(:,simInfo.bindingFlagsWrong) = new_rates(:,simInfo.bindingFlagsWrong)*simInfo.wrongFactorConcentration;
  end
    
  if equilibrium_flag
      % NL: for now I will implement a solution that assumes symmetry
      % between right and wrong cycles. Will need to revise if I choose to
      % examine more general architectures
      
      % generate input cell 
      forward_indices = simInfo.forwardRatesRight;
      backward_indices = simInfo.backwardRatesRight;
      
      % use rejection sampling to find reverse rates
      new_rates_temp = new_rates;
      err_flags = true(size(new_rates,1),1);
  
      % attempt to resample reverse rates         
      flux_factor = prod(new_rates_temp(err_flags,forward_indices),2) ./ prod(new_rates_temp(err_flags,backward_indices),2);
      flux_factor = flux_factor.^(1/length(backward_indices));
      new_rates_temp(err_flags,backward_indices) = new_rates_temp(err_flags,backward_indices).*flux_factor;          

      % find cases where this renormalization pushes rates outside of bounds
      err_flags = max(log10(new_rates_temp(:,backward_indices))<rate_bounds(1,backward_indices)...
        |log10(new_rates_temp(:,backward_indices))>rate_bounds(2,backward_indices),[],2);

      % try the reverse operation
      new_rates_temp(err_flags,backward_indices) = new_rates(err_flags,backward_indices); % reset
      new_rates_temp(err_flags,forward_indices) = new_rates_temp(err_flags,forward_indices)./flux_factor(err_flags);             

      err_flags = max(log10(new_rates_temp(:,forward_indices))<rate_bounds(1,forward_indices)...
        |log10(new_rates_temp(:,forward_indices))>rate_bounds(2,forward_indices),[],2);      
      
      % apply to wrong cycle
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

  % apply half max constraint if ncessary
  if simInfo.half_max_flag
      if true%~simInfo.numTesting
          if simInfo.nStates == 4
              new_rates = applyHMConstraintFourState(simInfo,new_rates,simInfo.rate_bounds);
          else
              new_rates = applyHMConstraintSixState(simInfo,new_rates,simInfo.rate_bounds);
          end
      else
          if simInfo.nStates == 4
              new_rates = applyHMConstraintFourStateNum(simInfo,new_rates,simInfo.rate_bounds);
          else
              error('Numerical HM not yet implemented for 6+ states')
          end
      end
              
          
  end
  

  
  
  
