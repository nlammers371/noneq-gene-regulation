
function [networkInfo, RSymSmall, RSym] = makeActivatorGeneralNetworkNS(networkInfo6,nSites,nGeneralReactions)
                             
    % extract key parameters
    n_right_bound_init = networkInfo6.n_right_bound;
    n_wrong_bound_init = networkInfo6.n_wrong_bound;
    activity_vec_init =networkInfo6.activeStateFilter==1;
    transitionInfoInit = networkInfo6.transitionInfoArray;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % call recursive loop to buil up transition rate array
    n_g_engaged = double(activity_vec_init);
    n_right_bound = n_right_bound_init;
    n_wrong_bound = n_wrong_bound_init;
    transitionInfoArray = transitionInfoInit;       
        
    % add cross-plane connections. Let us assume the 6 state plane that is
    % equivalent to the simpler 6 state model considered correspond to the
    % first block  
    
    % add binding reactions first
    for n = 1:nSites-1
      
        baseNew = length(n_g_engaged);
        rate_mask = repelem(eye(3),baseNew,baseNew);
        
        % expand transition array            
        transitionInfoArray = repmat(transitionInfoArray,3,3);
        transitionInfoArray = transitionInfoArray.*rate_mask;
        
        % simply triple ns engaged size every time
        n_g_engaged = repmat(n_g_engaged,1,3);

        % double and increment binding info vec
        n_right_bound = repmat(n_right_bound,1,3);   
        n_right_bound(baseNew+1:2*baseNew) = n_right_bound(baseNew+1:2*baseNew) + 1;
                 
        n_wrong_bound = repmat(n_wrong_bound,1,3);
        n_wrong_bound(2*baseNew+1:3*baseNew) = n_wrong_bound(2*baseNew+1:3*baseNew) + 1;
        
        for i = 1:baseNew
            binding_ids = [1 2];
            for j = 1:2
                ind1 = i;
                ind2 = baseNew*j + i;

                % add binding info
                transitionInfoArray(ind2,ind1) = binding_ids(j); 
                transitionInfoArray(ind1,ind2) = -binding_ids(j); 
            end
        end
    end
    
    % now general factor reactions
    for n = 1:nGeneralReactions-1
      
        baseNew = length(n_g_engaged);
        rate_mask = repelem(eye(2),baseNew,baseNew);
        
        % expand transition array            
        transitionInfoArray = repmat(transitionInfoArray,2,2);
        transitionInfoArray = transitionInfoArray.*rate_mask;
        
        % double and increment ns binding info
        n_g_engaged = repmat(n_g_engaged,1,2);  
        n_g_engaged(baseNew+1:end) = n_g_engaged(baseNew+1:end) + 1;       

        % double binding info vec
        n_right_bound = repmat(n_right_bound,1,2);   
        n_wrong_bound = repmat(n_wrong_bound,1,2);   
    
        for i = 1:baseNew

            ind1 = i;
            ind2 = baseNew + i;

            % add binding info
            transitionInfoArray(ind2,ind1) = 3; 
            transitionInfoArray(ind1,ind2) = -3; 

        end
    end
    
    % tallie the new stats
    n_total_bound = n_right_bound + n_wrong_bound;
    nStates = length(n_total_bound);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Generate symbolic transition rate matrix
   
    % initialize cooperativity factors
    new_sym_vec = sym({'wmp', 'wa', 'wi'},'positive');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % specify var order, var bounds, and whether it can be swept by default
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    networkInfo = struct;
    networkInfo.sweepVarList = [networkInfo6.sweepVarList new_sym_vec];
    networkInfo.sweepVarStrings = cellstr(string(networkInfo.sweepVarList));
    networkInfo.defaultValues = ones(size(networkInfo.sweepVarList));
    networkInfo.sweepFlags = ones(size(networkInfo.sweepVarList))==1;
    networkInfo.sweepFlags(1:3) = false;
    networkInfo.paramBounds = repmat([-5 ; 5], 1, length(networkInfo.sweepFlags));        
    networkInfo.cr_index = 1;
    networkInfo.cw_index = find(strcmp(networkInfo.sweepVarStrings,'cw'));
    networkInfo.a_index = find(strcmp(networkInfo.sweepVarStrings,'a'));
    
    % initialize symbolic variables
    for s = 1:length(networkInfo.sweepVarList)
        eval(['syms ' networkInfo.sweepVarStrings{s} ' positive;']);
    end
    
    % reduce bounds on cooperativity terms
%     coop_flags = contains(networkInfo.sweepVarStrings,'w')&~strcmp(networkInfo.sweepVarStrings,'cw');
%     networkInfo.paramBounds(2,coop_flags) = 3;
%     networkInfo.paramBounds(1,coop_flags) = -3;
    
    % adjust bounds to ensure activating behavior
    wip_index = contains(networkInfo.sweepVarStrings,'wip');
    networkInfo.paramBounds(2,wip_index) = 0;
    
    wap_index = contains(networkInfo.sweepVarStrings,'wap');
    networkInfo.paramBounds(1,wap_index) = 0;
    
    % realistic limits on binding rates
%     kp_index = contains(networkInfo.sweepVarStrings,'kp');
%     networkInfo.paramBounds(2,kp_index) = 2;
%     
%     wpa_index = contains(networkInfo.sweepVarStrings,'wpa');
%     networkInfo.paramBounds(2,wpa_index) = 0;
    
    % Provide info for equilibrium constraint application and cycle flux
    % calculations    
    networkInfo.forwardRateConstants(1) = {[wap wma]};
    networkInfo.forwardRateIndices(1) = {find(ismember(networkInfo.sweepVarList,networkInfo.forwardRateConstants{1}))};
    networkInfo.backwardRateConstants(1) = {[wip wpa]};
    networkInfo.backwardRateIndices(1) = {find(ismember(networkInfo.sweepVarList,networkInfo.backwardRateConstants{1}))}; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% perform basic assignments based on type
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    RSym = sym(zeros(nStates));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (1) add basic binding/unbinding reactions and NS reactions       
    binding_reactions = ismember(transitionInfoArray,[1,2]);
    unbinding_reactions = ismember(transitionInfoArray,-[1,2]);
    general_tf_engages = transitionInfoArray==3;
    general_tf_disengages = transitionInfoArray==-3;
    
    RSym(binding_reactions) = kp;
    RSym(unbinding_reactions) = km;
    
    RSym(general_tf_engages) = ka;
    RSym(general_tf_disengages) = ki;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (2) add NS->activator interactions
    for n = 1:nGeneralReactions
        RSym(binding_reactions&n_g_engaged==n) = RSym(binding_reactions&n_g_engaged==n)*wpa^n;        

        RSym(unbinding_reactions&n_g_engaged==n) = RSym(unbinding_reactions&n_g_engaged==n)*wma^n;        
    end   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (3) add activator->NS
    for n = 1:nSites
        RSym(general_tf_engages&n_total_bound==n) = RSym(general_tf_engages&n_total_bound==n)*wap^n;        

        RSym(general_tf_disengages&n_total_bound==n) = RSym(general_tf_disengages&n_total_bound==n)*wip^n;        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (4) add activator <-> activator   
    for n = 2:nSites
%          RSym(binding_reactions&n_total_bound==n-1) = RSym(binding_reactions&n_total_bound==n-1)*wpp^(n-1);
      
        RSym(unbinding_reactions&n_total_bound==n) = RSym(unbinding_reactions&n_total_bound==n)*wmp^(n-1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (5) add general <-> general   
    for n = 2:nGeneralReactions
        RSym(general_tf_engages&n_g_engaged==n-1) = RSym(general_tf_engages&n_g_engaged==n-1)*wa^(n-1);
        
        RSym(general_tf_disengages&n_g_engaged==n) = RSym(general_tf_disengages&n_g_engaged==n)*wi^(n-1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (6) add cr-dependence, cw-dependence, and affinity factor
    RSym(transitionInfoArray==1) = RSym(transitionInfoArray==1)*cr;
    RSym(transitionInfoArray==2) = RSym(transitionInfoArray==2)*cw;
    RSym(transitionInfoArray==-2) = RSym(transitionInfoArray==-2)*a;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (7) condense duplicate binding states
    full_state_array = [n_g_engaged' n_right_bound' n_wrong_bound'];
    [unique_state_array,ia,~] = unique(full_state_array,'rows');
    
    RSymSmall = sym(zeros(size(unique_state_array,1)));
    for i = 1:size(unique_state_array,1)        
        state_ids = ismember(full_state_array,unique_state_array(i,:),'rows');%n_g_engaged==unique_state_array(i,1) & n_right_bound==unique_state_array(i,2) & n_wrong_bound==unique_state_array(i,3);
        to_col = sum(RSym(state_ids,:),1);
        RSymSmall(i,:) = to_col(ia);
    end        
    
    % add diagonal factors 
    RSym(eye(size(RSym,1))==1) = -sum(RSym); 
    RSymSmall(eye(size(RSymSmall,1))==1) = -sum(RSymSmall); 
    
    % generate activity vectors
    n_g_engaged_small = n_g_engaged(ia);
    activity_vec_small = n_g_engaged_small==max(n_g_engaged_small);        
    activity_vec_full = n_g_engaged==max(n_g_engaged);       
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
    networkInfo.n_total_bound_full = n_total_bound;
    networkInfo.n_g_engaged_full = n_g_engaged;   
    networkInfo.n_right_bound = n_right_bound(ia);  
    networkInfo.n_wrong_bound = n_wrong_bound(ia);  
    networkInfo.n_wrong_bound_full = n_wrong_bound;  
    networkInfo.n_g_engaged = n_g_engaged_small;   
    networkInfo.n_total_bound = n_total_bound(ia);
    networkInfo.activeStateFilterFull = activity_vec_full;
    networkInfo.activeStateFilter = activity_vec_small;