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

    % generate flags indicating states where specific site is bound by
    % correct and incorrect factors    
    right_state_flags = repmat(n_right_bound_init,1,3)==1;
    wrong_state_flags = repmat(n_wrong_bound_init,1,3)==1;
    
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
    % Generate symbolic transition rate matrix

    % initialize core rate multiplier variables
    syms cr cw alphaa positive

    % initialize core locus activity transition rates
    syms kip kap kim kam positive

    % initialize core binding rates
    syms kpi kpa kmi kma positive

    % initialize weights to allow for impact of 2 bound
    syms kip2 kap2 kma2 kmi2 positive

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % specify var order, var bounds, and whether it can be swept by default
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    networkInfo = struct;
    networkInfo.sweepVarList = [cr cw alphaa kip kap kim kam kpi kmi kpa kma kap2 kip2 kma2 kmi2];
    networkInfo.sweepVarStrings = cellstr(string(networkInfo.sweepVarList));
    networkInfo.defaultValues = [1 100 100 1 1 1 1 1 1 1 1 1 1 1 1];
    networkInfo.sweepFlags = [0 0 0 1 1 1 1 1 1 1 1 1 1 1 1]==1;
    networkInfo.fourStateInputFlags = [1 0 0 1 1 1 1 1 1 1 1 0 0 0 0]==1;
    networkInfo.fourStateWrongInputFlags = [0 1 1 1 1 1 1 1 1 1 1 0 0 0 0]==1;
    networkInfo.bindingFlags = [0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 ];
    networkInfo.unbindingFlags = [0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 ];
    networkInfo.paramBounds = repmat([-4 ; 4 ], 1, length(networkInfo.sweepFlags));
%     networkInfo.paramBounds(:,end-3:end) = [-1 -1 -1 -1 ; 1 1 1 1];
    
    % require that activation/deactivation multiplers by greater/less than
    % 1. This is consistent with our assumptions that factors function as
    % activators
%     activationFlags = contains(networkInfo.sweepVarStrings,'wap');    
%     networkInfo.paramBounds(1,activationFlags) = 0;
%     deactivationFlags = contains(networkInfo.sweepVarStrings,'wip');
%     networkInfo.paramBounds(2,deactivationFlags) = 0;
    
    % zero out non-rate edges (NL: not sure this is strictly necessary)
    networkInfo.paramBounds(:,1:3) = 0;
    
    % networkInfo.paramBounds(1,end-3:end) = 0;
    networkInfo.cr_index = 1;
    networkInfo.cw_index = 2;
    networkInfo.b_index = 3;

    % Provide info for equilibrium constraint application and cycle flux
    % calculations

    % for 18 state network, we have to start with cooperativity and binding
    % factors
    networkInfo.forwardRateConstants(1) = {[kpi kma]};
    networkInfo.forwardRateAdjustFlags(1) = {[1 1]==1};
    networkInfo.forwardRateIndices(1) = {find(ismember(networkInfo.sweepVarList,networkInfo.forwardRateConstants{1}))};
    networkInfo.backwardRateConstants(1) = {[kpa kmi]};
    networkInfo.backwardRateAdjustFlags(1) = {[1 1]==1};
    networkInfo.backwardRateIndices(1) = {find(ismember(networkInfo.sweepVarList,networkInfo.backwardRateConstants{1}))};

    networkInfo.forwardRateConstants(2) = {[kap kim]};
    networkInfo.forwardRateAdjustFlags(2) = {[1 1]==1};
    networkInfo.forwardRateIndices(2) = {find(ismember(networkInfo.sweepVarList,networkInfo.forwardRateConstants{2}))};
    networkInfo.backwardRateConstants(2) = {[ kip kam]};
    networkInfo.backwardRateAdjustFlags(2) = {[1 1]==1};
    networkInfo.backwardRateIndices(2) = {find(ismember(networkInfo.sweepVarList,networkInfo.backwardRateConstants{2}))};

    networkInfo.forwardRateConstants(3) = {[kap2 kma2]};
    networkInfo.forwardRateAdjustFlags(3) = {[1 1]==1};
    networkInfo.forwardRateIndices(3) = {find(ismember(networkInfo.sweepVarList,networkInfo.forwardRateConstants{3}))};
    networkInfo.backwardRateConstants(3) = {[kip2 kmi2]};
    networkInfo.backwardRateAdjustFlags(3) = {[1 1]==1};
    networkInfo.backwardRateIndices(3) = {find(ismember(networkInfo.sweepVarList,networkInfo.backwardRateConstants{3}))};

    % HM constraints
    networkInfo.hmRatePairs = [kpi kpa; kmi kma; kip kim; kap kam]; 
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
    RSym(ismember(transitionInfoArray,[1,2]) & ~activity_vec_full) = kpi;
    RSym(ismember(transitionInfoArray,[1,2]) & activity_vec_full) = kpa;

    RSym(ismember(transitionInfoArray,-[1,2]) & ~activity_vec_full & n_total_bound<2) = kmi;
    RSym(ismember(transitionInfoArray,-[1,2]) & activity_vec_full & n_total_bound<2) = kma;
    
    RSym(ismember(transitionInfoArray,-[1,2]) & ~activity_vec_full & n_total_bound==2) = kmi2;
    RSym(ismember(transitionInfoArray,-[1,2]) & activity_vec_full & n_total_bound==2) = kma2;

    % basic locus fluctuation rates
    % assume that the number bound (but not the identity) can alter this rate
    RSym(ismember(transitionInfoArray,3) & n_total_bound==0) = kam;
    RSym(ismember(transitionInfoArray,3) & n_total_bound==1) = kap;
    RSym(ismember(transitionInfoArray,3) & n_total_bound==2) = kap2;

    RSym(ismember(transitionInfoArray,-3) & n_total_bound==0) = kim;
    RSym(ismember(transitionInfoArray,-3) & n_total_bound==1) = kip;
    RSym(ismember(transitionInfoArray,-3) & n_total_bound==2) = kip2;

    % layer on 2-bound multipliers for unbinding 
%     f1 = ismember(transitionInfoArray,-[1,2]) & n_total_bound==2 & ~activity_vec_full;
%     RSym(f1) = RSym(f1) * wmi2;
% 
%     f2 = ismember(transitionInfoArray,-[1,2]) & n_total_bound==2 & activity_vec_full;
%     RSym(f2) = RSym(f2) * wma2;

    % add specificity and concentration factors
    RSym(transitionInfoArray==1) = RSym(transitionInfoArray==1)*cr;
    RSym(transitionInfoArray==2) = RSym(transitionInfoArray==2)*cw;
    RSym(transitionInfoArray==-2) = RSym(transitionInfoArray==-2)*alphaa;

    % add diagonal factors 
    RSym(eye(size(RSym,1))==1) = -sum(RSym);

    % wrte to file
%     RSymFun = matlabFunction(RSym,'File',[savePath 'RSymFun'],'Optimize',true,...
%               'Vars',networkInfo.sweepVarList);                   
    
    % save helper variables
    activeStates = find(activity_vec_full);
    networkInfo.nStates = nStates;
    networkInfo.right_state_flags = right_state_flags;
    networkInfo.wrong_state_flags = wrong_state_flags;
    networkInfo.activeStates = activeStates;
    networkInfo.permittedConnections = permittedConnections;
    networkInfo.transitionInfoArray = transitionInfoArray;
    networkInfo.RSym = RSym;
    networkInfo.n_right_bound = n_right_bound;
    networkInfo.n_wrong_bound = n_wrong_bound;
    networkInfo.n_total_bound = n_total_bound;
    networkInfo.activeStateFilter = activity_vec_full;