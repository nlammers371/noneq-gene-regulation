function [networkInfo, RSym] = generate18StateNetwork_v2(baseNum, n_wrong_bound_init...
          ,n_right_bound_init,activity_vec_init,...
          transitionInfoInit,permittedConnectionsInit)

    nStates = 18;
    rate_mask = repelem(eye(3),baseNum,baseNum);

    % full connection array
    permittedConnections = repmat(permittedConnectionsInit,3,3);%permittedConnectionsInit
    permittedConnections = permittedConnections.*rate_mask;

    % full binding array
    transitionInfoArray = repmat(transitionInfoInit,3,3);%permittedConnectionsInit
    transitionInfoArray = transitionInfoArray.*rate_mask;

    % full activity vector
    activity_vec_full = repmat(activity_vec_init,1,3);

    % binding info vecs
    n_right_bound = repmat(n_right_bound_init,1,3);
    n_right_bound(baseNum+1:2*baseNum) = n_right_bound(baseNum+1:2*baseNum) + 1;

    n_wrong_bound = repmat(n_wrong_bound_init,1,3);
    n_wrong_bound(2*baseNum+1:3*baseNum) = n_wrong_bound(2*baseNum+1:3*baseNum) + 1;

    n_total_bound = n_wrong_bound + n_right_bound;
    
    %%
    % add cross-plane connections. Let us assume the 6 state plane that is
    % equivalent to the simpler 6 state model considered correspond to the
    % first block
    for i = 1:baseNum
        binding_ids = [1 2];
        for j = 1:2
            ind1 = i;
            ind2 = baseNum*j + i;

            % add to array
            permittedConnections(ind1,ind2) = 1;
            permittedConnections(ind2,ind1) = 1;

            % add binding info
            transitionInfoArray(ind2,ind1) = binding_ids(j); % specify whether it is a cognate or non-cognate factor binding
            transitionInfoArray(ind1,ind2) = -2; % all unbinding events from non-specific site are equivalent
        end
    end
    permittedConnections = permittedConnections==1;
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Generate symbolic transition rate matrix

    % initialize core rate multiplier variables
    syms cr cw alphaa positive

    % initialize core locus activity transition rates
    syms wip wap kim kam positive

    % initialize core binding rates
    syms kpi kpa kmi kma positive

    % initialize weights to allow for impact of 2 bound
    syms wma wmi positive

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % specify var order, var bounds, and whether it can be swept by default
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    networkInfo = struct;
    networkInfo.sweepVarList = [cr cw alphaa wip wap kim kam kpi kmi kpa kma wma wmi];
    networkInfo.sweepVarStrings = cellstr(string(networkInfo.sweepVarList));
    networkInfo.defaultValues = [1 1 100 1 1 1 1 1 1 1 1 1 1 1];
    networkInfo.sweepFlags = [0 0 0 1 1 1 1 1 1 1 1 1 1]==1;
    networkInfo.bindingFlags = [0 0 0 0 0 0 0 1 0 1 0 0 0];
    networkInfo.unbindingFlags = [0 0 0 0 0 0 0 0 1 0 1 0 0];
    networkInfo.paramBounds = repmat([-3 ; 2], 1, length(networkInfo.sweepFlags));
    
    % apply more stringent constraints to activator-NS interaction terms to
    % ensure activating TF activity
    wip_index = strcmp(networkInfo.sweepVarStrings,'wip');
    networkInfo.paramBounds(2,wip_index) = 0;
    
    wap_index = strcmp(networkInfo.sweepVarStrings,'wap');
    networkInfo.paramBounds(1,wap_index) = 0;
    
    % networkInfo.paramBounds(1,end-3:end) = 0;
    networkInfo.cr_index = 1;
    networkInfo.cw_index = 2;
    networkInfo.b_index = 3;
    
    % Provide info for equilibrium constraint application and cycle flux
    % calculations
    networkInfo.forwardRateConstants(1) = {[kpi kma wap]};
    networkInfo.forwardRateAdjustFlags(1) = {[1 1 1]==1};
    networkInfo.forwardRateIndices(1) = {find(ismember(networkInfo.sweepVarList,networkInfo.forwardRateConstants{1}))};
    networkInfo.backwardRateConstants(1) = {[kpa kmi wip]};
    networkInfo.backwardRateAdjustFlags(1) = {[1 1 1]==1};
    networkInfo.backwardRateIndices(1) = {find(ismember(networkInfo.sweepVarList,networkInfo.backwardRateConstants{1}))};

%     networkInfo.forwardRateConstants(2) = {[ kim]};
%     networkInfo.forwardRateAdjustFlags(2) = {[1 1]==1};
%     networkInfo.forwardRateIndices(2) = {find(ismember(networkInfo.sweepVarList,networkInfo.forwardRateConstants{2}))};
%     networkInfo.backwardRateConstants(2) = {[  kam]};
%     networkInfo.backwardRateAdjustFlags(2) = {[1 1]==1};
%     networkInfo.backwardRateIndices(2) = {find(ismember(networkInfo.sweepVarList,networkInfo.backwardRateConstants{2}))};

    networkInfo.forwardRateConstants(2) = {[wma]};
    networkInfo.forwardRateAdjustFlags(2) = {[1]==1};
    networkInfo.forwardRateIndices(2) = {find(ismember(networkInfo.sweepVarList,networkInfo.forwardRateConstants{2}))};
    networkInfo.backwardRateConstants(2) = {[wmi]};
    networkInfo.backwardRateAdjustFlags(2) = {[1]==1};
    networkInfo.backwardRateIndices(2) = {find(ismember(networkInfo.sweepVarList,networkInfo.backwardRateConstants{2}))};

    % HM constraints
    networkInfo.hmRatePairs = [kpi kpa; kmi kma; wip kim; wap kam]; 
    networkInfo.hmRatePairIndices = NaN(size(networkInfo.hmRatePairs));
    for i = 1:size(networkInfo.hmRatePairIndices,1)
        networkInfo.hmRatePairIndices(i,:) = find(ismember(networkInfo.sweepVarList,networkInfo.hmRatePairs(i,:)));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform basic assignments based on type
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    RSym = sym(zeros(nStates));

    % basic binding and unbinding rates (right and wrong factors have same base
    % rate)
    binding_flags = ismember(transitionInfoArray,[1,2]) ;
    RSym(binding_flags & ~activity_vec_full) = kpi;
    RSym(binding_flags & activity_vec_full) = kpa;

    unbinding_flags = ismember(transitionInfoArray,[-1,-2]);
    RSym(unbinding_flags & ~activity_vec_full & n_total_bound<2) = kmi;
    RSym(unbinding_flags & activity_vec_full & n_total_bound<2) = kma;
    
    RSym(unbinding_flags & ~activity_vec_full & n_total_bound==2) = kmi*wmi;
    RSym(unbinding_flags & activity_vec_full & n_total_bound==2) = kma*wma;

    % basic locus fluctuation rates
    % assume that the number bound (but not the identity) can alter this rate
    RSym(ismember(transitionInfoArray,3) & n_total_bound==0) = kam;
    RSym(ismember(transitionInfoArray,3) & n_total_bound==1) = kam*wap;
    RSym(ismember(transitionInfoArray,3) & n_total_bound==2) = kam*wap^2;

    RSym(ismember(transitionInfoArray,-3) & n_total_bound==0) = kim;
    RSym(ismember(transitionInfoArray,-3) & n_total_bound==1) = kim*wip;
    RSym(ismember(transitionInfoArray,-3) & n_total_bound==2) = kim*wip^2;

    % add specificity and concentration factors
    RSym(transitionInfoArray==1) = RSym(transitionInfoArray==1)*cr;
    RSym(transitionInfoArray==2) = RSym(transitionInfoArray==2)*cw;
    RSym(transitionInfoArray==-2) = RSym(transitionInfoArray==-2)*alphaa;

    % add diagonal factors 
    RSym(eye(size(RSym,1))==1) = -sum(RSym);  

    % save helper variables
    activeStates = find(activity_vec_full);
    networkInfo.nStates = nStates;
    networkInfo.RSym = RSym;
    networkInfo.activeStates = activeStates;
    networkInfo.permittedConnections = permittedConnections;
    networkInfo.transitionInfoArray = transitionInfoArray;
    networkInfo.n_right_bound = n_right_bound;    
    networkInfo.n_wrong_bound = n_wrong_bound;
    networkInfo.n_total_bound = n_total_bound;
    networkInfo.activeStateFilter = activity_vec_full;