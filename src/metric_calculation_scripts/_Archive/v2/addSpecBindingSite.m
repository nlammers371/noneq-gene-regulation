function [networkInfo, RSym] = addSpecBindingSite(networkInfoInit)

    % extract basic info to get started
    permittedConnectionsInit = networkInfoInit.permittedConnections;
    baseNum = size(permittedConnectionsInit,1);
    transitionInfoInit = networkInfoInit.transitionInfoArray;
    activity_vec_init = networkInfoInit.activeStateFilter;    
    n_right_bound_init = networkInfoInit.n_right_bound;
    
    nStatesOut = 2*baseNum; 
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
    
    % initialize weights to allow for impact of 2 bound
    for n = 3:nSites
        a = sym(['kap' num2str(nSites)],'positive');
        b = sym(['kip' num2str(nSites)],'positive');        
        c = sym(['kma' num2str(nSites)],'positive');    
        d = sym(['kmi' num2str(nSites)],'positive');
    end    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % specify var order, var bounds, and whether it can be swept by default
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    networkInfo = networkInfoInit;
    networkInfo.sweepVarList = [networkInfo.sweepVarList a b c d];
    networkInfo.sweepVarStrings = cellstr(string(networkInfo.sweepVarList));
    networkInfo.defaultValues = [networkInfo.defaultValues 1 1 1 1];
    networkInfo.sweepFlags = [networkInfo.sweepFlags true(1,4)];
    networkInfo.fourStateInputFlags = [networkInfo.fourStateInputFlags false(1,4)];
%     networkInfo.fourStateWrongInputFlags = [networkInfo.fourStateWrongInputFlags false(1,4)];
    networkInfo.bindingFlags = [0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 ];
    networkInfo.unbindingFlags = [0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 ];
    networkInfo.paramBounds = [networkInfo.paramBounds [-4 -4 -4 -4 ; 4 4 4 4]];
    
    % require that activation/deactivation multiplers by greater/less than
    % 1. This is consistent with our assumptions that factors function as
    % activators      

    % get index of previous state weights thta need to be included in new
    % normalization sets    
    fwd_prev_sym = sym(['kma' num2str(nSites-1)],'positive');
    bkd_prev_sym = sym(['kmi' num2str(nSites-1)],'positive');   
    
    % Add info for equilibrium constraint application and cycle flux
    % calculations    
    networkInfo.forwardRateConstants(end+1) = {[fwd_prev_sym a c]};
    networkInfo.forwardRateAdjustFlags(end+1) = {[0 1 1]==1};
    networkInfo.forwardRateIndices(end+1) = {find(ismember(networkInfo.sweepVarList,networkInfo.forwardRateConstants{end}))};
    networkInfo.backwardRateConstants(end+1) = {[bkd_prev_sym b d]};
    networkInfo.backwardRateAdjustFlags(end+1) = {[0 1 1]==1};
    networkInfo.backwardRateIndices(end+1) = {find(ismember(networkInfo.sweepVarList,networkInfo.backwardRateConstants{end}))};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate symbolic transition rate matrix

    % initialize core rate multiplier variables
    syms cr positive

    % initialize core locus activity transition rates
    syms kip kap kim kam positive

    % initialize core binding rates
    syms kpi kpa kmi kma positive

    % initialize weights to allow for impact of 2 bound
    syms kip2 kap2 kma2 kmi2 positive
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform basic assignments based on type
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    RSym = sym(zeros(nStatesOut));

    % basic binding and unbinding rates (right and wrong factors have same base
    % rate)
    RSym(transitionInfoArray==1 & ~activity_vec_full) = kpi;
    RSym(transitionInfoArray==1 & activity_vec_full) = kpa;

    RSym(transitionInfoArray==-1 & ~activity_vec_full & n_total_bound==1) = kmi;
    RSym(transitionInfoArray==-1 & ~activity_vec_full & n_total_bound==2) = kmi2;
    RSym(transitionInfoArray==-1 & activity_vec_full & n_total_bound==1) = kma;
    RSym(transitionInfoArray==-1 & activity_vec_full & n_total_bound==2) = kma2;

    % basic locus fluctuation rates
    % assume that the number bound (but not the identity) can alter this rate
    RSym(ismember(transitionInfoArray,3) & n_total_bound==0) = kam;
    RSym(ismember(transitionInfoArray,3) & n_total_bound==1) = kap;
    RSym(ismember(transitionInfoArray,3) & n_total_bound==2) = kap2;

    RSym(ismember(transitionInfoArray,-3) & n_total_bound==0) = kim;
    RSym(ismember(transitionInfoArray,-3) & n_total_bound==1) = kip;
    RSym(ismember(transitionInfoArray,-3) & n_total_bound==2) = kip2;

    % add higher order terms
    for n = 3:nSites
        % re-initialize symbolic variables
        a = sym(['kap' num2str(nSites)],'positive');
        b = sym(['kip' num2str(nSites)],'positive');        
        c = sym(['kma' num2str(nSites)],'positive');    
        d = sym(['kmi' num2str(nSites)],'positive');
        
        % activity treansition multipliers
        f1 = ismember(transitionInfoArray,3) & n_total_bound==n;
        RSym(f1) = a;
        
        f2 = ismember(transitionInfoArray,-3) & n_total_bound==n;
        RSym(f2) = b;
        
        % unbinding multipliers
        f3 = transitionInfoArray==-1 & n_total_bound>=n & ~activity_vec_full;
        RSym(f3) = d;

        f4 = transitionInfoArray==-1 & n_total_bound>=n & activity_vec_full;
        RSym(f4) = c;
    end     

    % add specificity and concentration factors
    RSym(transitionInfoArray==1) = RSym(transitionInfoArray==1)*cr;
    
    % add diagonal factors 
    RSym(eye(size(RSym,1))==1) = -sum(RSym);
    
    % save helper variables
    activeStates = find(activity_vec_full);
    networkInfo.nStates = nStatesOut;
    networkInfo.activeStates = activeStates;
    networkInfo.permittedConnections = permittedConnections;
    networkInfo.transitionInfoArray = transitionInfoArray;
    networkInfo.n_right_bound = n_right_bound;
    networkInfo.n_total_bound = n_total_bound;
    networkInfo.activeStateFilter = activity_vec_full;