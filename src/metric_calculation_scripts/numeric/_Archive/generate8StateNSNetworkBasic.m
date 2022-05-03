% this function generates key network infor for numeric simulations of
% networks with two non-spefic steps. For simplicity, we assume that each
% step has the same base parameters, but that ns factors can intereact
% with one another and with the transcription factor
function [networkInfo, RSym] = generate8StateNSNetworkBasic(networkInfo6)
                             
    % extract key parameters
    baseNum = 4;
    n_right_bound_init = networkInfo6.n_right_bound(1:4);
    activity_vec_init =networkInfo6.activeStateFilter(1:4)==1;
    transitionInfoInit = networkInfo6.transitionInfoArray(1:4,1:4);
                              
                              
    rate_mask = repelem(eye(2),baseNum,baseNum);
    nStates = size(rate_mask,1);

    % full binding array
    transitionInfoArray = repmat(transitionInfoInit,2,2);
    transitionInfoArray = transitionInfoArray.*rate_mask;

    % full activity vector
    activity_vec_full = [zeros(size(activity_vec_init)) activity_vec_init]==1;

    % binding info vecs
    n_right_bound = repmat(n_right_bound_init,1,2);      
    n_total_bound = n_right_bound;    
    
    %%
    % add cross-plane connections. Let us assume the 6 state plane that is
    % equivalent to the simpler 6 state model considered correspond to the
    % first block
    for i = 1:baseNum
     
        ind1 = i;
        ind2 = baseNum + i;

        % add binding info
        transitionInfoArray(ind2,ind1) = 4; % "4" signifies binding of a new non-specific factor
        transitionInfoArray(ind1,ind2) = -4; 
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Generate symbolic transition rate matrix

    % initialize core rate multiplier variables
    syms cr positive

    % initialize core locus activity transition rates
    syms wip wap ki ka positive
%     syms wip2 wap2 ki2 ka2 positive
    
    % initialize core unbinding rates
    syms km wma positive

    % binding rates
    syms kp wpa 
    
    % initialize weights to capture cooperativity between NS factors
    syms wa12 positive
    syms wi12 positive

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % specify var order, var bounds, and whether it can be swept by default
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    networkInfo = struct;
    networkInfo.sweepVarList = [cr wip wap ki ka km wma kp wpa wa12 wi12];
    networkInfo.sweepVarStrings = cellstr(string(networkInfo.sweepVarList));
    networkInfo.defaultValues = ones(size(networkInfo.sweepVarList));
    networkInfo.sweepFlags = ones(size(networkInfo.sweepVarList))==1;
    networkInfo.sweepFlags(1) = false;
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% perform basic assignments based on type
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    RSym = sym(zeros(nStates));

    % (1) add base reactions first
    RSym(transitionInfoArray==1) = kp;
    RSym(transitionInfoArray==-1) = km;
    
    RSym(transitionInfoArray==3) = ka;
    RSym(transitionInfoArray==-3) = ki;
    
    RSym(transitionInfoArray==4) = ka;
    RSym(transitionInfoArray==-4) = ki;
    
    % (2) layer on cooperative effects between activator and NS reactions
    for n = 3:4
        a = sym('wap','positive');
        b = sym('wip','positive');   
        c = sym('wma','positive');
        d = sym('wpa' ,'positive');
        
        % activator effects on NS
        ns_on_filter = transitionInfoArray==n&n_total_bound==1;
        ns_off_filter = transitionInfoArray==-n&n_total_bound==1;
        RSym(ns_on_filter) = RSym(ns_on_filter)*a;
        RSym(ns_off_filter) = RSym(ns_off_filter)*b;      
        
        % NS effects on activator 
        a_on_filter = transitionInfoArray==1&any(transitionInfoArray==-n);
        a_off_filter = transitionInfoArray==-1&any(transitionInfoArray==-n);
        RSym(a_on_filter) = RSym(a_on_filter)*d;
        RSym(a_off_filter) = RSym(a_off_filter)*c;              
    end
    
    % (3) Add NS-NS interactions
    
    % activator effects on NS
    on12_filter = transitionInfoArray==3&any(transitionInfoArray==-4);
    on21_filter = transitionInfoArray==4&any(transitionInfoArray==-3);
    RSym(on12_filter|on21_filter) = RSym(on12_filter|on21_filter)*wa12;
%     RSym(on21_filter) = RSym(on21_filter)*wa21; 
    
    off12_filter = transitionInfoArray==-3&any(transitionInfoArray==-4);
    off21_filter = transitionInfoArray==-4&any(transitionInfoArray==-3);
    RSym(off12_filter|off21_filter) = RSym(off12_filter|off21_filter)*wi12;
%     RSym(off21_filter) = RSym(off21_filter)*wi21; 
    
    % add concentration factors
    RSym(transitionInfoArray==1) = RSym(transitionInfoArray==1)*cr;

    % add diagonal factors 
    RSym(eye(size(RSym,1))==1) = -sum(RSym);  

    % save helper variables
    activeStates = find(activity_vec_full);
    networkInfo.nStates = nStates;
    networkInfo.RSym = RSym;
    networkInfo.activeStates = activeStates;
    networkInfo.transitionInfoArray = transitionInfoArray;
    networkInfo.n_right_bound = n_right_bound;    
    networkInfo.n_total_bound = n_total_bound;
    networkInfo.activeStateFilter = activity_vec_full;