% this function generates key network infor for numeric simulations of
% networks with two non-spefic steps. For simplicity, we assume that each
% step has the same base parameters, but that ns factors can intereact
% with one another and with the transcription factor
function [networkInfo, RSym] = addNSReaction(networkInfoInit)
                             
    % extract key parameters
    baseNum = size(networkInfoInit.RSym,1);
    RSymInit = networkInfoInit.RSym;
    n_right_bound_init = networkInfoInit.n_right_bound;
    activity_vec_init =networkInfoInit.activeStateFilter;
    transitionInfoInit = networkInfoInit.transitionInfoArray;
                              
    % Adding an NS reaction increases total number of states by 2
    rate_mask = repelem(eye(2),baseNum,baseNum);
    nStates = size(rate_mask,1);

    % full binding array
    transitionInfoArray = repmat(transitionInfoInit,2,2);
    transitionInfoArray = transitionInfoArray.*rate_mask;
    RSymRaw = repmat(RSymInit,2,2);
    RSymRaw = RSymRaw.*rate_mask;
    
    % full activity vector
    activity_vec_full = [zeros(size(activity_vec_init)) activity_vec_init]==1;

    % binding info vecs
    n_right_bound = repmat(n_right_bound_init,1,2);      
    n_total_bound = n_right_bound;    
    
    %%
    max_num = max(abs(transitionInfoInit(:)))+1;
    % add cross-plane connections
    for i = 1:baseNum
     
        ind1 = i;
        ind2 = baseNum + i;

        % add binding info
        transitionInfoArray(ind2,ind1) = max_num; 
        transitionInfoArray(ind1,ind2) = -max_num; 
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Generate symbolic transition rate matrix

    % infer number of NS reactions 
    n_ns = max_num-2;    
    sym_vec = sym(zeros(1,6));
    
    % initialize basal activation/inactivation variables
    sym_vec(1) = sym(['ka' num2str(n_ns)],'positive');       
    sym_vec(2) = sym(['ki' num2str(n_ns)],'positive'); 
    
    % initialize interactions NS->activator
    sym_vec(3) = sym(['wpa' num2str(n_ns)],'positive');
    sym_vec(4) = sym(['wma' num2str(n_ns)],'positive');
    
    % initialize interactions NS<-activator
    sym_vec(5) = sym(['wap' num2str(n_ns)],'positive');
    sym_vec(6) = sym(['wip' num2str(n_ns)],'positive');
        
    % initialize pairwise interaction variables between NS reactions   
    for i = 1:n_ns-1
%         ind = 4*i;
        % i->j
        sym_vec(end+1) = sym(['wa' num2str(n_ns) num2str(i)],'positive');
        sym_vec(end+1) = sym(['wi' num2str(n_ns) num2str(i)],'positive');
        
        % j -> i
        sym_vec(end+1) = sym(['wa' num2str(i) num2str(n_ns)],'positive');
        sym_vec(end+1) = sym(['wi' num2str(i) num2str(n_ns)],'positive');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % specify var order, var bounds, and whether it can be swept by default
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    networkInfo = struct;
    networkInfo.sweepVarList = [networkInfoInit.sweepVarList sym_vec];
    networkInfo.sweepVarStrings = cellstr(string(networkInfo.sweepVarList));
    networkInfo.defaultValues = ones(size(networkInfo.sweepVarList));
    networkInfo.sweepFlags = ones(size(networkInfo.sweepVarList))==1;
    networkInfo.sweepFlags(1) = false; % never sweep cr
    
    networkInfo.paramBounds = repmat([-6 ; 8], 1, length(networkInfo.sweepFlags));        
    networkInfo.cr_index = 1;
    
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
    % ...
    % NL: hold on this for now

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% perform basic assignments based on type
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    RSym = RSymRaw;%sym(zeros(nStates));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (1) add base reactions first    
    RSym(transitionInfoArray==max_num) = sym_vec(1);
    RSym(transitionInfoArray==-max_num) = sym_vec(2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (2) layer on cooperative effects between activator and NS     
    a = sym(['wap' num2str(n_ns)],'positive');
    b = sym(['wip' num2str(n_ns)],'positive');   
    c = sym(['wma' num2str(n_ns)],'positive');
    d = sym(['wpa' num2str(n_ns)],'positive');

    % activator -> NS
    ns_on_filter = transitionInfoArray==max_num&n_total_bound==1;
    ns_off_filter = transitionInfoArray==-max_num&n_total_bound==1;
    RSym(ns_on_filter) = RSym(ns_on_filter)*a;
    RSym(ns_off_filter) = RSym(ns_off_filter)*b;      

    % NS -> activator 
    a_on_filter = transitionInfoArray==1&any(transitionInfoArray==-max_num);
    a_off_filter = transitionInfoArray==-1&any(transitionInfoArray==-max_num);
    RSym(a_on_filter) = RSym(a_on_filter)*d;
    RSym(a_off_filter) = RSym(a_off_filter)*c;              
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (3) Add NS-NS interactions
    
    for n = 3:max_num-1
      
        % i->j
        a = sym(['wa' num2str(n_ns) num2str(n-2)],'positive');
        b = sym(['wi' num2str(n_ns) num2str(n-2)],'positive');
        
        % j -> i
        c = sym(['wa' num2str(n-2) num2str(n_ns)],'positive');
        d = sym(['wi' num2str(n-2) num2str(n_ns)],'positive');
        
        % activation effects 
        on12_filter = transitionInfoArray==n&any(transitionInfoArray==-max_num);
        on21_filter = transitionInfoArray==max_num&any(transitionInfoArray==-n);
        RSym(on12_filter) = RSym(on12_filter)*c;
        RSym(on21_filter) = RSym(on21_filter)*a; 
    
        off12_filter = transitionInfoArray==-n&any(transitionInfoArray==-max_num);
        off21_filter = transitionInfoArray==-max_num&any(transitionInfoArray==-n);
        RSym(off12_filter) = RSym(off12_filter)*d;
        RSym(off21_filter) = RSym(off21_filter)*b; 
    end
    
    % add concentration factors
%     RSym(transitionInfoArray==1) = RSym(transitionInfoArray==1)*cr;

    % add diagonal factors 
    RSym(eye(size(RSym,1))==1) = 0;  
    RSym(eye(size(RSym,1))==1) = -sum(RSym);  

    % save helper variables
    activeStates = find(activity_vec_full);
    networkInfo.nStates = nStates;
    networkInfo.RSym = RSym;
    networkInfo.activeStates = activeStates;
%     networkInfo.permittedConnections = permittedConnections;
    networkInfo.transitionInfoArray = transitionInfoArray;
    networkInfo.n_right_bound = n_right_bound;    
    networkInfo.n_total_bound = n_total_bound;
    networkInfo.activeStateFilter = activity_vec_full;