% this function generates key network infor for numeric simulations of
% networks with two non-spefic steps. For simplicity, we assume that each
% step has the sambe base parameters, but that ns factors can intereact
% with one another and with the transcription factor
function [networkInfo, RSym] = generate8StateNSDiffNetwork(baseNum...
                                ,n_right_bound_init,activity_vec_init,...
                                transitionInfoInit,permittedConnectionsInit)
    
    rate_mask = repelem(eye(2),baseNum,baseNum);
    nStates = size(rate_mask,1);
    
    % full connection array
    permittedConnections = repmat(permittedConnectionsInit,2,2);%permittedConnectionsInit
    permittedConnections = permittedConnections.*rate_mask;

    % full binding array
    transitionInfoArray = repmat(transitionInfoInit,2,2);%permittedConnectionsInit
    transitionInfoArray = transitionInfoArray.*rate_mask;

    % full activity vector
    activity_vec_full = [zeros(size(activity_vec_init)) activity_vec_init]==1;

    % binding info vecs
    n_right_bound = repmat(n_right_bound_init,1,2);      
    n_total_bound = n_right_bound;
    
    % binding infor for NS factors
    n_ns_bound_init = activity_vec_init*1;
    n_ns_bound = repmat(n_ns_bound_init,1,2);
    n_ns_bound(baseNum+1:end) = n_ns_bound(baseNum+1:end) + 1;
    
    %%
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
        transitionInfoArray(ind2,ind1) = 4; % "4" signifies binding of a new non-specific factor
        transitionInfoArray(ind1,ind2) = -4; % all unbinding events from non-specific site are equivalent
        
    end
    permittedConnections = permittedConnections==1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Generate symbolic transition rate matrix

    % initialize core rate multiplier variables
    syms cr positive

    % initialize core locus activity transition rates
    syms kip kap kim kam positive

    % initialize core binding rates
    syms kpi kpa kmi kma positive

    % initialize rates for NS interactions
    syms kip2 kap2 kim2 kam2 positive
    
    % initialize weights to capture cooperativity between factors
    for i = 3:4
    end
    syms kpa2 kma2 positive

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % specify var order, var bounds, and whether it can be swept by default
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    networkInfo = struct;
    networkInfo.sweepVarList = [cr kip kap kim kam kpi kmi kpa kma kip2 kap2 kim2 kam2 kpa2 kma2];
    networkInfo.sweepVarStrings = cellstr(string(networkInfo.sweepVarList));
    networkInfo.defaultValues = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
    networkInfo.sweepFlags = [0 1 1 1 1 1 1 1 1 1 1 1 1 1 1]==1;
    networkInfo.fourStateInputFlags = [1 1 1 1 1 1 1 1 1 0 0 0 0 0 0]==1;
    networkInfo.bindingFlags = [0 0 0 0 0 1 0 1 0 0 0 0 0 0 0];
    networkInfo.unbindingFlags = [0 0 0 0 0 0 1 0 1 0 0 0 0 0 0];
    networkInfo.paramBounds = repmat([-4 ; 4 ], 1, length(networkInfo.sweepFlags));
    
    % networkInfo.paramBounds(1,end-3:end) = 0;
    networkInfo.cr_index = 1;

    % Provide info for equilibrium constraint application and cycle flux
    % calculations
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

%     networkInfo.forwardRateConstants(3) = {[kap2 kma2]};
%     networkInfo.forwardRateAdjustFlags(3) = {[1 1]==1};
%     networkInfo.forwardRateIndices(3) = {find(ismember(networkInfo.sweepVarList,networkInfo.forwardRateConstants{3}))};
%     networkInfo.backwardRateConstants(3) = {[kip2 kmi2]};
%     networkInfo.backwardRateAdjustFlags(3) = {[1 1]==1};
%     networkInfo.backwardRateIndices(3) = {find(ismember(networkInfo.sweepVarList,networkInfo.backwardRateConstants{3}))};

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
    RSym(transitionInfoArray==1 & n_ns_bound==0) = kpi;
    RSym(transitionInfoArray==1 & n_ns_bound==1) = kpa;
    RSym(transitionInfoArray==1 & n_ns_bound==2) = kpa2;
    
    RSym(transitionInfoArray==-1 & n_ns_bound==0) = kmi;
    RSym(transitionInfoArray==-1 & n_ns_bound==1) = kma;
    RSym(transitionInfoArray==-1 & n_ns_bound==2) = kma2;   

    % basic locus fluctuation rates
    % assume that the number bound (but not the identity) can alter this rate
    RSym(ismember(transitionInfoArray,3) & n_total_bound==0) = kam;
    RSym(ismember(transitionInfoArray,3) & n_total_bound==1) = kap;
    
    RSym(ismember(transitionInfoArray,-3) & n_total_bound==0) = kim;
    RSym(ismember(transitionInfoArray,-3) & n_total_bound==1) = kip;

    RSym(ismember(transitionInfoArray,4) & n_total_bound==0) = kam2;
    RSym(ismember(transitionInfoArray,4) & n_total_bound==1) = kap2;
    
    RSym(ismember(transitionInfoArray,-4) & n_total_bound==0) = kim2;
    RSym(ismember(transitionInfoArray,-4) & n_total_bound==1) = kip2;
    
    % add concentration factors
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