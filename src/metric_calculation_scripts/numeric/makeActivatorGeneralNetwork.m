
function [networkInfo, RSymSmall, RSym] = makeActivatorGeneralNetwork(networkInfo6,nSpecificSites,nGeneralReactions)
                             
    % extract key parameters
    n_right_bound_init = networkInfo6.n_right_bound(1:4);
    activity_vec_init =networkInfo6.activeStateFilter(1:4)==1;
    transitionInfoInit = networkInfo6.transitionInfoArray(1:4,1:4);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % call recursive loop to buil up transition rate array
    n_g_engaged = double(activity_vec_init);
    n_right_bound = n_right_bound_init;
    transitionInfoArray = transitionInfoInit;
    
    % specific binding reactions first
    for n = 1:nSpecificSites-1
      
        baseNew = length(n_g_engaged);
        rate_mask = repelem(eye(2),baseNew,baseNew);
        
        % expand transition array            
        transitionInfoArray = repmat(transitionInfoArray,2,2);
        transitionInfoArray = transitionInfoArray.*rate_mask;
        
        % simply dobule ns engaged size every time
        n_g_engaged = repmat(n_g_engaged,1,2);

        % double and increment binding info vec
        n_right_bound = repmat(n_right_bound,1,2);      
        n_right_bound(baseNew+1:end) = n_right_bound(baseNew+1:end) + 1;            
    
        for i = 1:baseNew

            ind1 = i;
            ind2 = baseNew + i;

            % add binding info
            transitionInfoArray(ind2,ind1) = 1; 
            transitionInfoArray(ind1,ind2) = -1; 

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
    
        for i = 1:baseNew

            ind1 = i;
            ind2 = baseNew + i;

            % add binding info
            transitionInfoArray(ind2,ind1) = 3; 
            transitionInfoArray(ind1,ind2) = -3; 

        end
    end
    
    % tallie the new stats
    n_total_bound = n_right_bound;
    nStates = length(n_total_bound);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Generate symbolic transition rate matrix
    sym_vec = sym(zeros(1,0));
    % initialize cooperativity factors
    if nSpecificSites > 1
        sym_vec(1) = sym({'wmp'},'positive');
    end
    if nGeneralReactions > 1
        sym_vec(end+1:end+2) = sym({'wa', 'wi'},'positive');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % specify var order, var bounds, and whether it can be swept by default
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    networkInfo = struct;
    networkInfo.sweepVarList = [networkInfo6.sweepVarList([1 4:end]) sym_vec];
    networkInfo.sweepVarStrings = cellstr(string(networkInfo.sweepVarList));
    networkInfo.defaultValues = ones(size(networkInfo.sweepVarList));
    networkInfo.sweepFlags = ones(size(networkInfo.sweepVarList))==1;
    networkInfo.sweepFlags(1) = false;
    networkInfo.paramBounds = repmat([-5 ; 5], 1, length(networkInfo.sweepFlags));        
    networkInfo.cr_index = 1;
    
    % initialize symbolic variables
    for s = 1:length(networkInfo.sweepVarList)
        eval(['syms ' networkInfo.sweepVarStrings{s} ' positive;']);
    end
    
    % reduce bounds on cooperativity terms
%     coop_flags = contains(networkInfo.sweepVarStrings,'w');
%     networkInfo.paramBounds(2,coop_flags) = 3;
%     networkInfo.paramBounds(1,coop_flags) = -3;
    
    % adjust bounds to ensure activating behavior
    wip_index = contains(networkInfo.sweepVarStrings,'wip');
    networkInfo.paramBounds(2,wip_index) = 0;
    
    wap_index = contains(networkInfo.sweepVarStrings,'wap');
    networkInfo.paramBounds(1,wap_index) = 0;
    
    % impose realistic limits on binding rates
%     kp_index = contains(networkInfo.sweepVarStrings,'kp');
%     networkInfo.paramBounds(2,kp_index) = 2;
    
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
    RSym(transitionInfoArray==1) = kp;
    RSym(transitionInfoArray==-1) = km;
    
    RSym(transitionInfoArray==3) = ka;
    RSym(transitionInfoArray==-3) = ki;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (2) add NS->activator interactions
    for n = 1:nGeneralReactions
        RSym(transitionInfoArray==1&n_g_engaged==n) = RSym(transitionInfoArray==1&n_g_engaged==n)*wpa^n;        

        RSym(transitionInfoArray==-1&n_g_engaged==n) = RSym(transitionInfoArray==-1&n_g_engaged==n)*wma^n;        
    end   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (3) add activator->NS
    for n = 1:nSpecificSites
        RSym(transitionInfoArray==3&n_total_bound==n) = RSym(transitionInfoArray==3&n_total_bound==n)*wap^n;        

        RSym(transitionInfoArray==-3&n_total_bound==n) = RSym(transitionInfoArray==-3&n_total_bound==n)*wip^n;        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (4) add activator <-> activator   
    for n = 2:nSpecificSites
        RSym(transitionInfoArray==-1&n_total_bound==n) = RSym(transitionInfoArray==-1&n_total_bound==n)*wmp^(n-1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (5) add general <-> general   
    for n = 2:nGeneralReactions
        RSym(transitionInfoArray==3&n_g_engaged==n-1) = RSym(transitionInfoArray==3&n_g_engaged==n-1)*wa^(n-1);
        RSym(transitionInfoArray==-3&n_g_engaged==n) = RSym(transitionInfoArray==-3&n_g_engaged==n)*wi^(n-1);
    end           
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (6) add cr-dependence
    RSym(transitionInfoArray==1) = RSym(transitionInfoArray==1)*cr;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (7) condense duplicate binding states
    [unique_state_array,ia,~] = unique([n_g_engaged' n_total_bound'],'rows');
    
    RSymSmall = sym(zeros(size(unique_state_array,1)));
    for i = 1:size(unique_state_array,1)        
        state_ids = n_g_engaged==unique_state_array(i,1) & n_total_bound==unique_state_array(i,2);
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
    networkInfo.n_g_engaged = n_g_engaged_small;   
    networkInfo.n_total_bound = n_total_bound(ia);
    networkInfo.activeStateFilterFull = activity_vec_full;
    networkInfo.activeStateFilter = activity_vec_small;