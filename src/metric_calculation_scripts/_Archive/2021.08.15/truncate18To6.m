function networkInfo = truncate18To6(networkInfo18)

% filter key cvariable fields
networkInfo = networkInfo18;
keep_var_filter = ~contains(networkInfo18.sweepVarStrings,'wm');
networkInfo.sweepVarList = networkInfo18.sweepVarList(keep_var_filter);
networkInfo.sweepVarStrings = networkInfo18.sweepVarStrings(keep_var_filter);
networkInfo.defaultValues = networkInfo18.defaultValues(keep_var_filter);
networkInfo.sweepFlags = networkInfo18.sweepFlags(keep_var_filter);
networkInfo.bindingFlags = networkInfo18.bindingFlags(keep_var_filter);
networkInfo.unbindingFlags = networkInfo18.unbindingFlags(keep_var_filter);
networkInfo.paramBounds = networkInfo18.paramBounds(:,keep_var_filter);

% filter the network arrays
networkInfo.permittedConnections = networkInfo18.permittedConnections(1:6,1:6);
networkInfo.transitionInfoArray = networkInfo18.transitionInfoArray(1:6,1:6);

% adjust equilibrium constraints
networkInfo.forwardRateConstants = networkInfo18.forwardRateConstants(1);
networkInfo.forwardRateAdjustFlags = networkInfo18.forwardRateAdjustFlags(1);
networkInfo.forwardRateIndices = networkInfo18.forwardRateIndices(1);

networkInfo.backwardRateConstants = networkInfo18.backwardRateConstants(1);
networkInfo.backwardRateAdjustFlags = networkInfo18.backwardRateAdjustFlags(1);
networkInfo.backwardRateIndices = networkInfo18.backwardRateIndices(1);

% normalize RSym
RSym = networkInfo.RSym(1:6,1:6);
RSym(eye(size(RSym,1))==1) = 0;
RSym(eye(size(RSym,1))==1) = -sum(RSym);
networkInfo.RSym = RSym;

% other vector fields
networkInfo.n_total_bound = networkInfo18.n_total_bound(1:6);
networkInfo.n_right_bound = networkInfo18.n_right_bound(1:6);
networkInfo.n_wrong_bound = networkInfo18.n_wrong_bound(1:6);
networkInfo.activeStateFilter = networkInfo18.activeStateFilter(1:6)==1;