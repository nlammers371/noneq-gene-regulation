function new_rates = sample_rates(lb_array,ub_array,orig_rate_array,prop_sigma,...
          constrained_rates, half_max_flag, rate_bounds, cycleTime, c_val)
  
  
  if isempty(orig_rate_array)
    new_rates = reshape(10.^(prop_sigma*trandn(lb_array,ub_array)),[],8); % function to sample from truncated normal dist
  else
    new_rates = 10.^(orig_rate_array+reshape(prop_sigma*trandn(lb_array,ub_array),[],8));
  end   
       
  
  if ~isempty(constrained_rates)
    
    % helper variables
    k_filter = 1:8 < 5;
    r_filter = 1:8 > 4;
    flux_denom = sum(constrained_rates&k_filter);
    new_rates_temp = new_rates;
    
    % if necessary, renormalize rates to preserve detailed balance
    net_flux_vec = prod(new_rates(:,constrained_rates&k_filter),2) ./ prod(new_rates(:,constrained_rates&r_filter),2);    
    new_rates_temp(:,constrained_rates&r_filter) = net_flux_vec.^(1/flux_denom) .* new_rates(:,constrained_rates&r_filter);

    % find cases where this renormalization pushes rates outside of bounds
    err_flags1 = max(log10(new_rates_temp(:,constrained_rates&r_filter))<rate_bounds(1,constrained_rates&r_filter)...
      |log10(new_rates_temp(:,constrained_rates&r_filter))>rate_bounds(2,constrained_rates&r_filter),[],2);

    % reset
    new_rates_temp(err_flags1,:) = new_rates(err_flags1,:);
    new_rates_temp(err_flags1,constrained_rates&k_filter) = ...
      new_rates(err_flags1,constrained_rates&k_filter)./(net_flux_vec(err_flags1).^(1/flux_denom));
    
    err_flags2 = max(log10(new_rates_temp(:,constrained_rates&k_filter))<rate_bounds(1,constrained_rates&k_filter)...
      |log10(new_rates_temp(:,constrained_rates&k_filter))>rate_bounds(2,constrained_rates&k_filter),[],2);
    
    % if that didn't work, just set to NaN
    new_rates_temp(err_flags2,:) = NaN;
    new_rates = new_rates_temp;
    
  end
  
  if half_max_flag
    [new_rates,~] = constrainHalfMax(new_rates,1,rate_bounds);    
  elseif ~isempty(cycleTime)
    warning('NL: not sure this is implemented properly')
    new_rates = constrainCycleTime(new_rates, c_val, ...
                                         rate_bounds, cycleTime);
  end
  
  
  
