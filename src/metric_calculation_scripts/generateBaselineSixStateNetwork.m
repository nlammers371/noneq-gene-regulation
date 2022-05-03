function networkInfo = generateBaselineSixStateNetwork

    % define some basic parameters
    activeStatesFull = [3 4 5];
    baseNum = 6;

    % zero out forbidden transitions
    [fromRef,toRef] = meshgrid(1:baseNum,1:baseNum);
    diffRef = abs(fromRef-toRef);
    toRefHalved = toRef<=baseNum/2;
    fromRefHalved = fromRef<=baseNum/2;
    permittedConnectionsRaw= (diffRef==1 & toRefHalved==fromRefHalved) | diffRef==baseNum/2;

    % permute these connections to follow a more intuitive labeling scheme
    indexCurr = 1:baseNum;
    indexAdjusted = circshift(indexCurr,floor(baseNum/4));
    indexAdjusted = [indexAdjusted(1:baseNum/2) fliplr(indexAdjusted(baseNum/2+1:end))];
    [~,si] = sort(indexAdjusted);
    permittedConnections = permittedConnectionsRaw(si,si);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate an array with binding info
    transitionInfoArray = zeros(size(permittedConnections));

    % specific binding/unbinding
    spec_pairs = {[1,2],[4,3]};
    for s = 1:length(spec_pairs)
        ind1 = spec_pairs{s}(1);
        ind2 = spec_pairs{s}(2);
        % update
        transitionInfoArray(ind2,ind1) = 1;
        transitionInfoArray(ind1,ind2) = -1;
    end

    % non-specific binding/unbinding
    non_spec_pairs = {[1,6],[4,5]};
    for s = 1:length(non_spec_pairs)
        ind1 = non_spec_pairs{s}(1);
        ind2 = non_spec_pairs{s}(2);
        % update
        transitionInfoArray(ind2,ind1) = 2;
        transitionInfoArray(ind1,ind2) = -2;
    end

    % locus activity fluctuations
    locus_pairs = {[1,4],[6,5],[2,3]};
    for s = 1:length(locus_pairs)
        ind1 = locus_pairs{s}(1);
        ind2 = locus_pairs{s}(2);
        % update
        transitionInfoArray(ind2,ind1) = 3;
        transitionInfoArray(ind1,ind2) = -3;
    end

    % generate array tht indicates activity state
    activity_vec_full = false(1,baseNum);
    activity_vec_full(activeStatesFull) = 1;

    % generate flags indicating numbers of right and wrong factors bound
    n_right_bound = zeros(size(activity_vec_full));
    n_wrong_bound = n_right_bound;
    n_right_bound([2 3]) = 1;
    n_wrong_bound([5 6]) = 1;

    permittedConnections = permittedConnections==1;
    n_total_bound = n_wrong_bound + n_right_bound;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate symbolic transition rate matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % initialize core rate multiplier variables
    syms cr cw a positive

    % initialize baseline reaction rates
    syms ki ka km kp positive

    % initialize interaction terms
    syms wap wip wpa wma positive

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % specify var order, var bounds, and whether it can be swept by default
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    networkInfo = struct;
    networkInfo.sweepVarList = [cr cw a ki ka kp km wap wip wma wpa];
    networkInfo.sweepVarStrings = cellstr(string(networkInfo.sweepVarList));
    networkInfo.defaultValues = [1 1 100 1 1 1 1 1 1 1 1];
    networkInfo.sweepFlags = [0 0 0 1 1 1 1 1 1 1 1]==1;      
    networkInfo.paramBounds = repmat([-4 ; 4 ], 1, length(networkInfo.sweepFlags));    
    networkInfo.cr_index = 1;
    networkInfo.cw_index = 2;
    networkInfo.a_index = 3;

    % add constraints to ensure activating behavior
    wip_index = strcmp(networkInfo.sweepVarStrings,'wip');
    networkInfo.paramBounds(2,wip_index) = 0;
    
    wap_index = strcmp(networkInfo.sweepVarStrings,'wap');
    networkInfo.paramBounds(1,wap_index) = 0;
    
    % Provide info for equilibrium constraint application and cycle flux
    % calculations
    networkInfo.forwardRateConstants = [wap wma];
    networkInfo.forwardRateIndices = {find(ismember(networkInfo.sweepVarList,networkInfo.forwardRateConstants))};
    networkInfo.backwardRateConstants = [wip wpa];
    networkInfo.backwardRateIndices = {find(ismember(networkInfo.sweepVarList,networkInfo.backwardRateConstants))};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform basic assignments based on type
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    RSym = sym(zeros(baseNum));

    % basic binding and unbinding rates
    RSym(ismember(transitionInfoArray,[1,2]) & ~activity_vec_full) = kp;
    RSym(ismember(transitionInfoArray,[1,2]) & activity_vec_full) = kp*wpa;

    RSym(ismember(transitionInfoArray,-[1,2]) & ~activity_vec_full) = km;
    RSym(ismember(transitionInfoArray,-[1,2]) & activity_vec_full) = km*wma;

    % basic locus fluctuation rates
    RSym(ismember(transitionInfoArray,3) & n_total_bound==0) = ka;
    RSym(ismember(transitionInfoArray,3) & n_total_bound==1) = ka*wap;

    RSym(ismember(transitionInfoArray,-3) & n_total_bound==0) = ki;
    RSym(ismember(transitionInfoArray,-3) & n_total_bound==1) = ki*wip;

    % add specificity and concentration factors
    RSym(transitionInfoArray==1) = RSym(transitionInfoArray==1)*cr;
    RSym(transitionInfoArray==2) = RSym(transitionInfoArray==2)*cw;
    RSym(transitionInfoArray==-2) = RSym(transitionInfoArray==-2)*a;

    RSym(eye(size(RSym,1))==1) = -sum(RSym);

    % add fields to network inf
    networkInfo.RSym = RSym;
    activeStates = find(activity_vec_full);
    networkInfo.nStates = baseNum;
    networkInfo.RSym = RSym;
    networkInfo.activeStates = activeStates;
    networkInfo.permittedConnections = permittedConnections;
    networkInfo.transitionInfoArray = transitionInfoArray;
    networkInfo.n_right_bound = n_right_bound;    
    networkInfo.n_wrong_bound = n_wrong_bound;
    networkInfo.n_total_bound = n_total_bound;
    networkInfo.activeStateFilter = activity_vec_full==1;