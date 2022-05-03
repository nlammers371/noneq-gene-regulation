function [networkInfo, RSym] = addSpecBindingSite_v3(networkInfoInit)

    % extract basic info to get started
    permittedConnectionsInit = networkInfoInit.permittedConnections;
    baseNum = size(permittedConnectionsInit,1);
    transitionInfoInit = networkInfoInit.transitionInfoArray;
    activity_vec_init = networkInfoInit.activeStateFilter;    
    n_right_bound_init = networkInfoInit.n_right_bound;
    
    nStates = 2*baseNum; 
    rate_mask = repelem(eye(2),baseNum,baseNum);

    % full connection array
    permittedConnections = repmat(permittedConnectionsInit,2,2);%permittedConnectionsInit
    permittedConnections = permittedConnections.*rate_mask;

    % full binding array
    transitionInfoArray = repmat(transitionInfoInit,2,2);%permittedConnectionsInit
    transitionInfoArray = transitionInfoArray.*rate_mask;

    % full activity vector
    activity_vec_full = repmat(activity_vec_init,1,2);

    % binding info vecs
    n_right_bound = repmat(n_right_bound_init,1,2);
    n_right_bound(baseNum+1:2*baseNum) = n_right_bound(baseNum+1:2*baseNum) + 1;

    n_total_bound = n_right_bound;
  
    % add cross-plane connections. Let us assume the 6 state plane that is
    % equivalent to the simpler 6 state model considered correspond to the
    % first block
    for i = 1:baseNum
       
        ind1 = i;
        ind2 = baseNum + i;

        % add to array
        permittedConnections(ind1,ind2) = 1;
        permittedConnections(ind2,ind1) = 1;

        % add binding info
        transitionInfoArray(ind2,ind1) = 1; % specify whether it is a cognate or non-cognate factor binding
        transitionInfoArray(ind1,ind2) = -1; % all unbinding events from non-specific site are equivalent
      
    end
    permittedConnections = permittedConnections==1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate number of sites in final network
    nSites = max(n_total_bound);        

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % specify var order, var bounds, and whether it can be swept by default
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    networkInfo = networkInfoInit;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate symbolic transition rate matrix

    % initialize core rate multiplier variables
    syms cr positive

    % initialize core locus activity transition rates
    syms wip wap kim kam positive

    % initialize core binding rates
    syms kpi kpa kmi kma positive

    % initialize weights to allow for impact of 2 bound
    syms wma wmi positive    
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform basic assignments based on type
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    RSym = sym(zeros(nStates));

    % basic binding and unbinding rates (right and wrong factors have same base
    % rate)
    RSym(transitionInfoArray==1 & ~activity_vec_full) = kpi;
    RSym(transitionInfoArray==1 & activity_vec_full) = kpa;

    RSym(transitionInfoArray==-1 & ~activity_vec_full) = kmi;
    RSym(transitionInfoArray==-1 & activity_vec_full) = kma;
    
    % basic locus fluctuation rates
    % assume that the number bound (but not the identity) can alter this rate
    for n = 0:nSites        
        RSym(ismember(transitionInfoArray,3) & n_total_bound==n) = kam*wap^n;
        RSym(ismember(transitionInfoArray,-3) & n_total_bound==n) = kim*wip^n;
    end
    
    % layer on 2-bound multipliers for unbinding 
    for n = 2:nSites     
        wf = n*(n-1)/2;
        ub_flags = ismember(transitionInfoArray,-1) & n_total_bound==n;
        RSym(ub_flags & activity_vec_full) = RSym(ub_flags & activity_vec_full)*wf*wma;
        RSym(ub_flags & ~activity_vec_full) = RSym(ub_flags & ~activity_vec_full)*wf*wmi;
    end

    % add specificity and concentration factors
    RSym(transitionInfoArray==1) = RSym(transitionInfoArray==1)*cr;
    
    % add diagonal factors 
    RSym(eye(size(RSym,1))==1) = -sum(RSym);
    
    % save helper variables
    activeStates = find(activity_vec_full);
    networkInfo.nStates = nStates;
    networkInfo.activeStates = activeStates;
    networkInfo.permittedConnections = permittedConnections;
    networkInfo.transitionInfoArray = transitionInfoArray;
    networkInfo.n_right_bound = n_right_bound;
    networkInfo.n_total_bound = n_total_bound;
    networkInfo.activeStateFilter = activity_vec_full;