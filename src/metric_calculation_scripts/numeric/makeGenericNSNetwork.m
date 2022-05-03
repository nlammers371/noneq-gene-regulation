
function [networkInfo, RSymSmall, RSym] = makeGenericNSNetwork(networkInfo6,nReactions)
                             
    % extract key parameters
    n_right_bound_init = networkInfo6.n_right_bound(1:4);
    activity_vec_init =networkInfo6.activeStateFilter(1:4)==1;
    transitionInfoInit = networkInfo6.transitionInfoArray(1:4,1:4);
    
    % call recursive loop to buil up transition rate array
%     activity_vec_full = activity_vec_init;
    n_ns_engaged = double(activity_vec_init);
    n_right_bound = n_right_bound_init;
    transitionInfoArray = transitionInfoInit;
    
    for n = 1:nReactions-1
      
        baseNew = length(n_right_bound);
        rate_mask = repelem(eye(2),baseNew,baseNew);
        
        % expand transition array            
        transitionInfoArray = repmat(transitionInfoArray,2,2);
        transitionInfoArray = transitionInfoArray.*rate_mask;
        
        % full activity vector
%         activity_vec_full = repmat(activity_vec_full,1,2);

        % binding info vec
        n_right_bound = repmat(n_right_bound,1,2);  
        
        % ns info vec
        n_ns_engaged = repmat(n_ns_engaged,1,2);  
        n_ns_engaged(baseNew+1:end) = n_ns_engaged(baseNew+1:end) + 1;            
    
        for i = 1:baseNew

            ind1 = i;
            ind2 = baseNew + i;

            % add binding info
            transitionInfoArray(ind2,ind1) = 3; 
            transitionInfoArray(ind1,ind2) = -3; 

        end
    end
    n_activator_bound = n_right_bound;
    nStates = length(n_activator_bound);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Generate symbolic transition rate matrix
   
    % initialize cooperativity factor
    syms wa wi positive
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % specify var order, var bounds, and whether it can be swept by default
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    networkInfo = struct;
    networkInfo.sweepVarList = [networkInfo6.sweepVarList([1 4:end]) wa wi];
    networkInfo.sweepVarStrings = cellstr(string(networkInfo.sweepVarList));
    networkInfo.defaultValues = ones(size(networkInfo.sweepVarList));
    networkInfo.sweepFlags = ones(size(networkInfo.sweepVarList))==1;
    networkInfo.sweepFlags(1) = false;
    networkInfo.paramBounds = repmat([-6 ; 6], 1, length(networkInfo.sweepFlags));        
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
    
    for n = 1:nReactions
        RSym(transitionInfoArray==1&n_ns_engaged==n) = RSym(transitionInfoArray==1&n_ns_engaged==n)*wpa^n;        

        RSym(transitionInfoArray==-1&n_ns_engaged==n) = RSym(transitionInfoArray==-1&n_ns_engaged==n)*wma^n;        
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (3) add activator->NS    
    RSym(transitionInfoArray==3&n_activator_bound==1) = RSym(transitionInfoArray==3&n_activator_bound==1)*wap;        
    RSym(transitionInfoArray==-3&n_activator_bound==1) = RSym(transitionInfoArray==-3&n_activator_bound==1)*wip;        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (4) add ns <-> ns   
    for n = 2:nReactions
        RSym(transitionInfoArray==-3&n_ns_engaged==n) = RSym(transitionInfoArray==-3&n_ns_engaged==n)*wi^(n-1);
        RSym(transitionInfoArray==3&n_ns_engaged==n) = RSym(transitionInfoArray==3&n_ns_engaged==n)*wa^(n-1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (5) add cr-dependence
    RSym(transitionInfoArray==1) = RSym(transitionInfoArray==1)*cr;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (6) condense duplicate binding states
    ns_flags = any(transitionInfoArray==-3);
    [unique_state_array,ia,~] = unique([1*n_ns_engaged' n_activator_bound'],'rows');
    
    % the key is to sum all flows INTO a composite state
    RSymSmall = sym(zeros(size(unique_state_array,1)));
    for i = 1:size(unique_state_array,1)        
        state_ids = n_ns_engaged==unique_state_array(i,1) & n_activator_bound==unique_state_array(i,2);
        to_row = sum(RSym(state_ids,:),1);
        RSymSmall(i,:) = to_row(ia);
    end
    
    % add diagonal factors 
    RSym(eye(size(RSym,1))==1) = -sum(RSym); 
    RSymSmall(eye(size(RSymSmall,1))==1) = -sum(RSymSmall); 
    
    % generate activity vectors
    n_ns_engaged_small = n_ns_engaged(ia);
    activity_vec_small = n_ns_engaged_small==max(n_ns_engaged_small);        
    activity_vec_full = n_ns_engaged==max(n_ns_engaged);       
    activeStatesFull = find(activity_vec_full);
    activeStatesSmall = find(activity_vec_small);
    
    % save helper variables
    networkInfo.nStates = nStates;
    networkInfo.RSymFull = RSym;
    networkInfo.RSym = RSymSmall;
    networkInfo.activeStatesFull = activeStatesFull;
    networkInfo.activeStatesSmall = activeStatesSmall;
    networkInfo.transitionInfoArrayFull = transitionInfoArray;
    networkInfo.transitionInfoArray = transitionInfoArray(ia,ia);
    networkInfo.n_right_bound_full = n_right_bound;    
    networkInfo.n_total_bound_full = n_activator_bound;
    networkInfo.n_right_bound = n_right_bound(ia);    
    networkInfo.n_total_bound = n_activator_bound(ia);
    networkInfo.activeStateFilterFull = activity_vec_full;
    networkInfo.activeStateFilter = activity_vec_small;