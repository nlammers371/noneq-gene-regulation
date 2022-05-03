function new_params_sweep = sample_rates_multi_v2(lb_array,ub_array,orig_param_array,prop_sigma,simInfo)

  % extract key parameters
  paramBounds = simInfo.paramBounds;  
  equilibrium_flag = simInfo.equilibrium_flag;
  sweepFlags = simInfo.sweepFlags;
  
  % generate proposed rates
  if isempty(orig_param_array)
      new_params_sweep = reshape(10.^(prop_sigma*trandn(lb_array,ub_array)),[],sum(sweepFlags)); % function to sample from truncated normal dist
  else
      new_params_sweep = 10.^(orig_param_array+reshape(prop_sigma*trandn(lb_array,ub_array),[],sum(sweepFlags)));
  end   
    
  % set param values for static variables (need this for HM calculations)
  new_params = NaN(size(new_params_sweep,1),length(sweepFlags));
  new_params(:,sweepFlags) = new_params_sweep;
  new_params(:,~sweepFlags) = repmat(simInfo.defaultValues(~sweepFlags),size(new_params,1),1); 
  
  % deal with collapse options
%   new_params(:,simInfo.map_from_indices) = new_params(:,simInfo.map_to_indices);
  
  % apply equality constraints
  constrained_indices = [simInfo.constraintPairCell{:}];
  for c = 1:length(simInfo.constraintPairCell)
      c_indices = simInfo.constraintPairCell{c};
      new_params(:,c_indices(2:end)) = new_params(:,c_indices(1));            
  end

  if equilibrium_flag
      % NL: for now I will implement a solution that assumes symmetry
      % between right and wrong cycles. Will need to revise if I choose to
      % examine more general architectures
      
      % first check to see if there are flux factors, in which case we
      % simply must set these to 0
      if isfield(simInfo,'fluxFactorFlags')
          new_params(:,simInfo.fluxFactorFlags==1) = 1;
      else
          forward_index_cell= simInfo.forwardRateIndices;
          backward_index_cell = simInfo.backwardRateIndices;
          use_flags = true;
          if isfield(simInfo,'forwardRateAdjustFlags')
              forward_flag_cell = simInfo.forwardRateAdjustFlags;
              backward_flag_cell = simInfo.backwardRateAdjustFlags;          
          else
              use_flags = false;
          end
          % use rejection sampling to find reverse rates
          new_params_temp = new_params;
          err_flags = false(size(new_params,1),1);

          for i = 1:length(forward_index_cell)

              update_flags = ~err_flags;

              % generate input cell 
              forward_indices = forward_index_cell{i};
              backward_indices = backward_index_cell{i};         

              % first check whether any indices are to constrained variables
              forward_flags = ~ismember(forward_indices,constrained_indices);
              backward_flags = ~ismember(backward_indices,constrained_indices);

              % extract adjust flags
              if use_flags
                  forward_flags = forward_flag_cell{i} & forward_flags;
                  backward_flags = backward_flag_cell{i} & backward_flags;
              end

              % attempt to resample reverse rates         
              flux_factor = prod(new_params_temp(:,forward_indices),2) ./ prod(new_params_temp(:,backward_indices),2);
              flux_factor = flux_factor.^(1/sum(backward_flags));
              new_params_temp(update_flags,backward_indices(backward_flags)) = ...
                    new_params_temp(update_flags,backward_indices(backward_flags)).*flux_factor(update_flags);          

              % find cases where this renormalization pushes rates outside of bounds
              oob_flags = max(log10(new_params_temp(:,backward_indices))<paramBounds(1,backward_indices)...
                |log10(new_params_temp(:,backward_indices))>paramBounds(2,backward_indices),[],2);  

              % identify subset of out-of-bounds errors that pertain to
              % relevant networks
              retry_flags = oob_flags & update_flags;

              % try the reverse operation
              new_params_temp(retry_flags,backward_indices(backward_flags)) = new_params(retry_flags,backward_indices(backward_flags)); % reset
              new_params_temp(retry_flags,forward_indices(forward_flags)) = new_params_temp(retry_flags,forward_indices(forward_flags))./flux_factor(retry_flags);                                   

              oob_flags2 = max(log10(new_params_temp(:,forward_indices))<paramBounds(1,forward_indices)...
                |log10(new_params_temp(:,forward_indices))>paramBounds(2,forward_indices),[],2);

              err_flags = oob_flags2 & update_flags;

              % set all networks that cannot be made to fit bounds to NaN
              new_params_temp(err_flags,:) = NaN;
          end
          new_params = new_params_temp;  
      end
  end

  
  % apply half max constraint if ncessary
  if simInfo.half_max_flag
      if simInfo.numCalcFlag
          new_params = applyHMConstraintNumeric(new_params, simInfo);          
      else
          if ismember(simInfo.nStates, [4,6])
              new_params = applyHMConstraintSymbolic(new_params,simInfo);       
          else
              error('Symbolic HM constraint not currently supported for >6 states')
          end
      end
  end
  

  new_params_sweep = new_params(:,sweepFlags);
  
  
