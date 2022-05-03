
function [networkInfo, RSymSmall, RSym] = generate8StateActivatorNetwork(networkInfo6)
                             
    % extract key parameters
    baseNum = 4;
    n_right_bound_init = networkInfo6.n_right_bound(1:4);
    activity_vec_init =networkInfo6.activeStateFilter(1:4)==1;
    transitionInfoInit = networkInfo6.transitionInfoArray(1:4,1:4);
    
    RSymRaw = networkInfo6.RSym(1:4,1:4);
    RSymRaw(eye(size(RSymRaw,1))==1) = 0;
    
    rate_mask = repelem(eye(2),baseNum,baseNum);
    nStates = size(rate_mask,1);

    % full binding array
    transitionInfoArray = repmat(transitionInfoInit,2,2);
    transitionInfoArray = transitionInfoArray.*rate_mask;

    RSymRaw = repmat(RSymRaw,2,2);
    RSymRaw = RSymRaw.*rate_mask;
    
    % full activity vector
    activity_vec_full = repmat(activity_vec_init,1,2);

    % binding info vecs
    n_right_bound = repmat(n_right_bound_init,1,2);      
    n_right_bound(baseNum+1:end) = n_right_bound(baseNum+1:end) + 1;
    n_total_bound = n_right_bound;    
    
    for i = 1:baseNum
     
        ind1 = i;
        ind2 = baseNum + i;

        % add binding info
        transitionInfoArray(ind2,ind1) = 1; 
        transitionInfoArray(ind1,ind2) = -1; 
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Generate symbolic transition rate matrix
   
    % initialize cooperativity factor
    syms wmp positive
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % specify var order, var bounds, and whether it can be swept by default
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    networkInfo = struct;
    networkInfo.sweepVarList = [networkInfo6.sweepVarList([1 4:end]) wmp];
    networkInfo.sweepVarStrings = cellstr(string(networkInfo.sweepVarList));
    networkInfo.defaultValues = ones(size(networkInfo.sweepVarList));
    networkInfo.sweepFlags = ones(size(networkInfo.sweepVarList))==1;
    networkInfo.sweepFlags(1) = false;
    networkInfo.paramBounds = repmat([-6 ; 4], 1, length(networkInfo.sweepFlags));        
    networkInfo.cr_index = 1;
    
    % initialize symbolic variables
    for s = 1:length(networkInfo.sweepVarList)
        eval(['syms ' networkInfo.sweepVarStrings{s} ' positive;']);
    end
    
    % reduce bounds on cooperativity terms
    coop_flags = contains(networkInfo.sweepVarStrings,'w');
    networkInfo.paramBounds(2,coop_flags) = 3;
    networkInfo.paramBounds(1,coop_flags) = -3;
    
    % adjust bounds to ensure activating behavior
    wip_index = contains(networkInfo.sweepVarStrings,'wip');
    networkInfo.paramBounds(2,wip_index) = 0;
    
    wap_index = contains(networkInfo.sweepVarStrings,'wap');
    networkInfo.paramBounds(1,wap_index) = 0;
    
    % Provide info for equilibrium constraint application and cycle flux
    % calculations

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% perform basic assignments based on type
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    RSym = sym(zeros(nStates));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (1) add basic binding/unbinding reactions and NS reactions   
    RSym(transitionInfoArray==1) = kp;
    RSym(transitionInfoArray==-1) = km;
    
    RSym(transitionInfoArray==3) = ka;
    RSym(transitionInfoArray==-3) = ki;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (2) add NS->activator interactions
    ns_engaged_flags = any(transitionInfoArray==-3);
    RSym(transitionInfoArray==1&ns_engaged_flags) = RSym(transitionInfoArray==1&ns_engaged_flags)*wpa;
    RSym(transitionInfoArray==-1&ns_engaged_flags) = RSym(transitionInfoArray==-1&ns_engaged_flags)*wma;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (3) add activator->NS
    RSym(transitionInfoArray==3&n_total_bound==1) = RSym(transitionInfoArray==3&n_total_bound==1)*wap;
    RSym(transitionInfoArray==3&n_total_bound==2) = RSym(transitionInfoArray==3&n_total_bound==2)*wap^2;
    
    RSym(transitionInfoArray==-3&n_total_bound==1) = RSym(transitionInfoArray==-3&n_total_bound==1)*wip;
    RSym(transitionInfoArray==-3&n_total_bound==2) = RSym(transitionInfoArray==-3&n_total_bound==2)*wip^2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (4) add activator<->activator        
    RSym(transitionInfoArray==-1&n_total_bound==2) = RSym(transitionInfoArray==-1&n_total_bound==2)*wmp;
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (5) condense duplicate binding states
    ns_flags = any(transitionInfoArray==-3);
    [unique_state_array,ia,~] = unique([1*ns_flags' n_total_bound'],'rows');
    
    RSymSmall = sym(zeros(size(unique_state_array,1)));
    for i = 1:size(unique_state_array,1)        
        state_ids = ns_flags==unique_state_array(i,1) & n_total_bound==unique_state_array(i,2);
        to_col = sum(RSym(state_ids,:),1);
        RSymSmall(i,:) = to_col(ia);
    end
    
    % add diagonal factors 
    RSym(eye(size(RSym,1))==1) = -sum(RSym); 
    RSymSmall(eye(size(RSymSmall,1))==1) = -sum(RSymSmall); 
    activity_vec_small = activity_vec_full(ia);
    
    % save helper variables
    activeStatesFull = find(activity_vec_full);
    activeStatesSmall = find(activity_vec_small);
    networkInfo.nStates = nStates;
    networkInfo.RSymFull = RSym;
    networkInfo.RSym = RSymSmall;
    networkInfo.activeStatesFull = activeStatesFull;
    networkInfo.activeStatesSmall = activeStatesSmall;
%     networkInfo.permittedConnections = permittedConnections;
    networkInfo.transitionInfoArray = transitionInfoArray;
    networkInfo.n_right_bound = n_right_bound;    
    networkInfo.n_total_bound = n_total_bound;
    networkInfo.activeStateFilterFull = activity_vec_full;
    networkInfo.activeStateFilter = activity_vec_small;