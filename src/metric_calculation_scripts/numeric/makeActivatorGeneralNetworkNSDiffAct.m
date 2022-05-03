
function [networkInfo, RSymSmall, RSym] = makeActivatorGeneralNetworkNSDiffAct(networkInfo6,nSites,nGeneralReactions)
                             
    % extract key parameters
    n_right_bound_init = networkInfo6.n_right_bound;
    n_wrong_bound_init = networkInfo6.n_wrong_bound;
    activity_vec_init = networkInfo6.activeStateFilter==1;
    transitionInfoInit = networkInfo6.transitionInfoArray;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % call recursive loop to buil up transition rate array
    n_g_engaged = double(activity_vec_init);
    n_right_bound = n_right_bound_init;
    n_wrong_bound = n_wrong_bound_init;
%     n_total_bound = n_right_bound+n_wrong_bound;
    transitionInfoArray = transitionInfoInit;       
        
    % initialize ID arrays to track status of specific 

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
        
        % simply dobule ns engaged size every time
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
                transitionInfoArray(ind2,ind1) = binding_ids(j)+n/10; 
                transitionInfoArray(ind1,ind2) = -binding_ids(j)-n/10; 
            end
        end
    end
    
    % now general factor reactions
    % note that we allow for distinct characteristics for the general
    % factors
%     g_factor = 4;
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
            transitionInfoArray(ind2,ind1) = 3 + n/10; 
            transitionInfoArray(ind1,ind2) = -3 - n/10; 

        end
%         g_factor = g_factor + 1;
    end
    
    % tallie the new stats
    n_total_bound = n_right_bound + n_wrong_bound;
    nStates = length(n_total_bound);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Initialize symbols and generate transition rate matrix
    % initialize basic variables
    syms cr cw a kp km positive
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform basic assignments based on type
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    RSym = sym(zeros(nStates));

    % generate additional reference arrays
    transitionInfoArrayRound = floor(abs(transitionInfoArray)).*sign(transitionInfoArray);
    binding_reactions = ismember(transitionInfoArrayRound,[1,2]);
    unbinding_reactions = ismember(transitionInfoArrayRound,-[1,2]);
    
    % generate book-keeping array to keep track of which general reactions
    % are in play
    general_id_array = zeros(nGeneralReactions,length(n_total_bound));
    for n = 1:nGeneralReactions
        general_id_array(n,:) = any(transitionInfoArray==-(3+(n-1)/10));
    end        
    
    diff_array = transitionInfoArray-transitionInfoArrayRound;
    diff_array(~(unbinding_reactions|binding_reactions)) = 0;
    diff_array(unbinding_reactions) = round(diff_array(unbinding_reactions) - 0.1,2);
    diff_array(binding_reactions) = round(diff_array(binding_reactions) + 0.1,2);
    
    activator_id_array = zeros(nSites,length(n_total_bound));
    for n = 1:nSites
        a_mask = round(diff_array,2)==-(n/10);
        a_id = min((1*a_mask).*transitionInfoArrayRound);
        activator_id_array(n,:) = abs(a_id);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (1) add basic binding/unbinding reactions and NS reactions           
    
    RSym(binding_reactions) = kp;
    RSym(unbinding_reactions) = km;
   
    % Initialize vector to store symbolic variables
    sym_vec = sym(zeros(1,0));
        
    sym_vec(end+1) = sym('ka','positive');       
    sym_vec(end+1) = sym('ki','positive');  

    % add to RSym
    general_tf_engages = transitionInfoArrayRound==3;
    general_tf_disengages = transitionInfoArrayRound==-3;
    RSym(general_tf_engages) = sym_vec(end-1);
    RSym(general_tf_disengages) = sym_vec(end);

    
    % initialize pairwise interaction variables between Activators 
    for i = 1:nSites
        for j = i+1:nSites
            % i->j
%             sym_vec(end+1) = sym(['wa' num2str(j) num2str(i)],'positive');
            sym_vec(end+1) = sym(['wmp' num2str(j) num2str(i)],'positive');
            % j -> i
%             sym_vec(end+1) = sym(['wa' num2str(i) num2str(j)],'positive');
            sym_vec(end+1) = sym(['wmp' num2str(i) num2str(j)],'positive');
            
            % add to RSym
            %i->j
%             RSym(diff_array==(i-1)/10 & general_id_array(j,:)) = RSym(transitionInfoArray==3+(i-1)/10 & general_id_array(j,:))*sym_vec(end-3);
            RSym(diff_array==-i/10 & activator_id_array(j,:)>0) = RSym(diff_array==-i/10 & activator_id_array(j,:)>0)*sym_vec(end);
            %j->i
%             RSym(transitionInfoArray==3+(j-1)/10 & general_id_array(i,:)) = RSym(transitionInfoArray==3+(j-1)/10 & general_id_array(i,:))*sym_vec(end-1);
            RSym(diff_array==-j/10 & activator_id_array(i,:)>0) = RSym(diff_array==-j/10 & activator_id_array(i,:)>0)*sym_vec(end-1);
        end
    end
    
    % initialize pairwise interaction variables betweenGeneral factors
    % interactions
    syms wa wi positive
    sym_vec(end+1) = wa;
    sym_vec(end+1) = wi;
    for n = 2:nGeneralReactions
        RSym(general_tf_engages&n_g_engaged==n-1) = RSym(general_tf_engages&n_g_engaged==n-1)*wa^(n-1);
        
        RSym(general_tf_disengages&n_g_engaged==n) = RSym(general_tf_disengages&n_g_engaged==n)*wi^(n-1);
    end
          
    % activator <-> GS
    
    for j = 1:nSites

        % initialize interactions NS->activator
        sym_vec(end+1) = sym(['wpa' num2str(j)],'positive');
        sym_vec(end+1) = sym(['wma' num2str(j)],'positive');

        % initialize interactions NS<-activator
        sym_vec(end+1) = sym(['wap' num2str(j)],'positive');
        sym_vec(end+1) = sym(['wip' num2str(j)],'positive');

        % add to RSym
        % a->g
        RSym(general_tf_engages & activator_id_array(j,:)>0) = RSym(general_tf_engages & activator_id_array(j,:)>0)*sym_vec(end-1);
        RSym(general_tf_disengages & activator_id_array(j,:)>0) = RSym(general_tf_disengages & activator_id_array(j,:)>0)*sym_vec(end);

        % g->a
        for n = 1:nGeneralReactions             
            RSym(diff_array==j/10&n_g_engaged==n) = RSym(diff_array==j/10&n_g_engaged==n)*sym_vec(end-3)^n;
            RSym(diff_array==-j/10&n_g_engaged==n) = RSym(diff_array==-j/10&n_g_engaged==n)*sym_vec(end-2)^n;
        end        

    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (6) add cr-dependence, cw-dependence, and affinity factor        
    RSym(transitionInfoArrayRound==1) = RSym(transitionInfoArrayRound==1)*cr;
    RSym(transitionInfoArrayRound==2) = RSym(transitionInfoArrayRound==2)*cw;
    RSym(transitionInfoArrayRound==-2) = RSym(transitionInfoArrayRound==-2)*a;
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % specify var order, var bounds, and whether it can be swept by default
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    networkInfo = struct;
    networkInfo.sweepVarList = [cr cw a km kp sym_vec];
    networkInfo.sweepVarStrings = cellstr(string(networkInfo.sweepVarList));
    networkInfo.defaultValues = ones(size(networkInfo.sweepVarList));
    networkInfo.sweepFlags = ones(size(networkInfo.sweepVarList))==1;
    networkInfo.sweepFlags(1:3) = false;
    networkInfo.paramBounds = repmat([-8 ; 4], 1, length(networkInfo.sweepFlags));        
    networkInfo.cr_index = 1;
    networkInfo.cw_index = find(strcmp(networkInfo.sweepVarStrings,'cw'));
    networkInfo.a_index = find(strcmp(networkInfo.sweepVarStrings,'a'));
    
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
    
    % impose realistic limits on binding rates
    kp_index = contains(networkInfo.sweepVarStrings,'kp');
    networkInfo.paramBounds(2,kp_index) = 2;
    
    wpa_index = contains(networkInfo.sweepVarStrings,'wpa');
    networkInfo.paramBounds(2,wpa_index) = 0;
    
    % Provide info for equilibrium constraint application and cycle flux
    % calculations

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (7) condense duplicate binding states
    full_state_array = [n_g_engaged' activator_id_array'];
    [unique_state_array,ia,ic] = unique(full_state_array,'rows');
%     ia = 1:size(RSym,1);
    
    RSymSmall = sym(zeros(size(unique_state_array,1)));
    for i = 1:size(unique_state_array,1)        
        state_ids = ismember(full_state_array,unique_state_array(i,:),'rows');
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
    
   
    
    %% save helper variables
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
    networkInfo.g_factor_id_array_full = general_id_array;
    networkInfo.g_factor_id_array = general_id_array(:,ia);
    networkInfo.activeStateFilterFull = activity_vec_full;
    networkInfo.activeStateFilter = activity_vec_small;