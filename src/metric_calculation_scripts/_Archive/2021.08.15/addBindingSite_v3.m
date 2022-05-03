function [networkInfo, RSym] = addBindingSite_v3(networkInfoInit)

    % extract basic info to get started
    permittedConnectionsInit = networkInfoInit.permittedConnections;
    baseNum = size(permittedConnectionsInit,1);
    transitionInfoInit = networkInfoInit.transitionInfoArray;
    activity_vec_init = networkInfoInit.activeStateFilter;
    n_wrong_bound_init = networkInfoInit.n_wrong_bound;
    n_right_bound_init = networkInfoInit.n_right_bound;
    
    nStatesOut = 3*baseNum;
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate number of sites in final network
    nSites = max(n_total_bound);
    
    % initialize weights to allow for impact of 2 bound
%     for n = 3:nSites
%         a = sym(['wap' num2str(nSites)],'positive');
%         b = sym(['wip' num2str(nSites)],'positive');        
%         c = sym(['wma' num2str(nSites)],'positive');    
%         d = sym(['wmi' num2str(nSites)],'positive');
%     end    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % specify var order, var bounds, and whether it can be swept by default
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    networkInfo = networkInfoInit;
%     networkInfo.sweepVarList = [networkInfo.sweepVarList a b c d];
%     networkInfo.sweepVarStrings = cellstr(string(networkInfo.sweepVarList));
%     networkInfo.defaultValues = [networkInfo.defaultValues 1 1 1 1];
%     networkInfo.sweepFlags = [networkInfo.sweepFlags true(1,4)];
%     networkInfo.fourStateInputFlags = [networkInfo.fourStateInputFlags false(1,4)];
%     networkInfo.fourStateWrongInputFlags = [networkInfo.fourStateWrongInputFlags false(1,4)];
%     networkInfo.bindingFlags = [0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 ];
%     networkInfo.unbindingFlags = [0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 ];
%     networkInfo.paramBounds = [networkInfo.paramBounds [-1 -1 -1 -1 ; 1 1 1 1]];
    
    % require that activation/deactivation multiplers by greater/less than
    % 1. This is consistent with our assumptions that factors function as
    % activators
%     activationFlags = contains(networkInfo.sweepVarStrings,'wap');    
%     networkInfo.paramBounds(1,activationFlags) = 0;
%     deactivationFlags = contains(networkInfo.sweepVarStrings,'wip');
%     networkInfo.paramBounds(2,deactivationFlags) = 0;            

    % get index of previous state weights that need to be included in new
    % normalization sets    
%     fwd_prev_sym = sym(['wma' num2str(nSites-1)],'positive');
%     bkd_prev_sym = sym(['wmi' num2str(nSites-1)],'positive');   
    
    % Add info for equilibrium constraint application and cycle flux
    % calculations    
%     networkInfo.forwardRateConstants(end+1) = {[fwd_prev_sym a c]};
%     networkInfo.forwardRateAdjustFlags(end+1) = {[0 1 1]==1};
%     networkInfo.forwardRateIndices(end+1) = {find(ismember(networkInfo.sweepVarList,networkInfo.forwardRateConstants{end}))};
%     networkInfo.backwardRateConstants(end+1) = {[bkd_prev_sym b d]};
%     networkInfo.backwardRateAdjustFlags(end+1) = {[0 1 1]==1};
%     networkInfo.backwardRateIndices(end+1) = {find(ismember(networkInfo.sweepVarList,networkInfo.backwardRateConstants{end}))};

%     % HM constraints
%     networkInfo.hmRatePairs = [kpi kpa; kmi kma; kip kim; kap kam]; 
%     networkInfo.hmRatePairIndices = NaN(size(networkInfo.hmRatePairs));
%     for i = 1:size(networkInfo.hmRatePairIndices,1)
%         networkInfo.hmRatePairIndices(i,:) = find(ismember(networkInfo.sweepVarList,networkInfo.hmRatePairs(i,:)));
%     end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate symbolic transition rate matrix

    % initialize core rate multiplier variables
    syms cr cw alphaa positive

    % initialize core locus activity transition rates
    syms kim kam positive

    % initialize core binding rates
    syms kpi kpa kmi kma positive

    % initialize weights to allow for impact of 2 bound
    syms wip wap wma wmi positive
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform basic assignments based on type
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    RSym = sym(zeros(nStatesOut));

    % basic binding and unbinding rates (right and wrong factors have same base
    % rate)
    RSym(ismember(transitionInfoArray,[1,2]) & ~activity_vec_full) = kpi;
    RSym(ismember(transitionInfoArray,[1,2]) & activity_vec_full) = kpa;

    RSym(ismember(transitionInfoArray,-[1,2]) & ~activity_vec_full) = kmi;
    RSym(ismember(transitionInfoArray,-[1,2]) & activity_vec_full) = kma;

    % basic locus fluctuation rates
    % assume that the number bound (but not the identity) can alter this rate
    for n = 0:nSites        
        RSym(ismember(transitionInfoArray,3) & n_total_bound==n) = kam*wap^n;
        RSym(ismember(transitionInfoArray,-3) & n_total_bound==n) = kim*wip^n;
    end
    % layer on 2-bound multipliers for unbinding 
    for n = 2:nSites     
        wf = n*(n-1)/2;
        ub_flags = ismember(transitionInfoArray,[-1,-2]) & n_total_bound==n;
        RSym(ub_flags & activity_vec_full) = kma*wf*wma;
        RSym(ub_flags & ~activity_vec_full) = kmi*wf*wmi;
    end
%     f1 = ismember(transitionInfoArray,-[1,2]) & n_total_bound>1 & ~activity_vec_full;
%     RSym(f1) = RSym(f1) * wmi2;
% 
%     f2 = ismember(transitionInfoArray,-[1,2]) & n_total_bound>1 & activity_vec_full;
%     RSym(f2) = RSym(f2) * wma2;

    % add specificity and concentration factors
    RSym(transitionInfoArray==1) = RSym(transitionInfoArray==1)*cr;
    RSym(transitionInfoArray==2) = RSym(transitionInfoArray==2)*cw;
    RSym(transitionInfoArray==-2) = RSym(transitionInfoArray==-2)*alphaa;

    % add higher order terms
%     for n = 3:nSites
%         % re-initialize symbolic variables
%         a = sym(['wap' num2str(nSites)],'positive');
%         b = sym(['wip' num2str(nSites)],'positive');        
%         c = sym(['wma' num2str(nSites)],'positive');    
%         d = sym(['wmi' num2str(nSites)],'positive');
%         
%         % activity treansition multipliers
%         f1 = ismember(transitionInfoArray,3) & n_total_bound>=n;
%         RSym(f1) = RSym(f1)*a;
%         
%         f2 = ismember(transitionInfoArray,-3) & n_total_bound>=n;
%         RSym(f2) = RSym(f2)*b;
%         
%         % unbinding multipliers
%         f3 = ismember(transitionInfoArray,-[1,2]) & n_total_bound>=n & ~activity_vec_full;
%         RSym(f3) = RSym(f3) * d;
% 
%         f4 = ismember(transitionInfoArray,-[1,2]) & n_total_bound>=n & activity_vec_full;
%         RSym(f4) = RSym(f4) * c;
%     end     

    % add diagonal factors 
    RSym(eye(size(RSym,1))==1) = -sum(RSym);
    
    % save helper variables
    activeStates = find(activity_vec_full);
    networkInfo.nStates = nStatesOut;
    networkInfo.activeStates = activeStates;
    networkInfo.RSym = RSym;
    networkInfo.permittedConnections = permittedConnections;
    networkInfo.transitionInfoArray = transitionInfoArray;
    networkInfo.n_right_bound = n_right_bound;
    networkInfo.n_wrong_bound = n_wrong_bound;
    networkInfo.n_total_bound = n_total_bound;
    networkInfo.activeStateFilter = activity_vec_full;