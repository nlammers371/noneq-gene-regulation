function networkInfo = truncate8To4(networkInfo8)

% filter key cvariable fields
networkInfo = networkInfo8;
keep_var_filter = ~contains(networkInfo8.sweepVarStrings,'wm');
networkInfo.sweepVarList = networkInfo8.sweepVarList(keep_var_filter);
networkInfo.sweepVarStrings = networkInfo8.sweepVarStrings(keep_var_filter);
networkInfo.defaultValues = networkInfo8.defaultValues(keep_var_filter);
networkInfo.sweepFlags = networkInfo8.sweepFlags(keep_var_filter);
networkInfo.bindingFlags = networkInfo8.bindingFlags(keep_var_filter);
networkInfo.unbindingFlags = networkInfo8.unbindingFlags(keep_var_filter);
networkInfo.paramBounds = networkInfo8.paramBounds(:,keep_var_filter);

% filter the network arrays
networkInfo.permittedConnections = networkInfo8.permittedConnections(1:4,1:4);
networkInfo.transitionInfoArray = networkInfo8.transitionInfoArray(1:4,1:4);

% adjust equilibrium constraints
networkInfo.forwardRateConstants = networkInfo8.forwardRateConstants(1);
networkInfo.forwardRateAdjustFlags = networkInfo8.forwardRateAdjustFlags(1);
networkInfo.forwardRateIndices = networkInfo8.forwardRateIndices(1);

networkInfo.backwardRateConstants = networkInfo8.backwardRateConstants(1);
networkInfo.backwardRateAdjustFlags = networkInfo8.backwardRateAdjustFlags(1);
networkInfo.backwardRateIndices = networkInfo8.backwardRateIndices(1);

% normalize RSym
RSym = networkInfo.RSym(1:4,1:4);
RSym(eye(size(RSym,1))==1) = 0;
RSym(eye(size(RSym,1))==1) = -sum(RSym);
networkInfo.RSym = RSym;

% other vector fields
networkInfo.n_total_bound = networkInfo8.n_total_bound(1:4);
networkInfo.n_right_bound = networkInfo8.n_right_bound(1:4);
networkInfo.activeStateFilter = networkInfo8.activeStateFilter(1:4)==1;