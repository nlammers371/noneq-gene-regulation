function [networkInfo, RSym] = generate6StateNucleosomeNetwork(n_wrong_bound...
          ,n_right_bound,activity_vec_full, transitionInfoArray, savePath)

    nStates = length(n_wrong_bound);        
    n_total_bound = n_wrong_bound + n_right_bound;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate symbolic transition rate matrix    
        
    % initialize core rate multiplier variables
    syms cr cw a positive % extra "a" on alpha is just to keep matlab 

    % initialize nucleosome on/off rates
    syms ki ka positive

    % initialize core binding rates
    syms kp km positive

    % initialize cooperativity weights % note that wi*wa = wp
    syms wi wa positive
    wp = wi*wa;
    % initialize non-eq factor (taken to be greater than 1)
    syms g positive
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % specify var order, var bounds, and whether it can be swept by default
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    networkInfo = struct;
    networkInfo.sweepVarList = [cr cw a ki ka kp km wi wa g];
    networkInfo.sweepVarStrings = cellstr(string(networkInfo.sweepVarList));
    networkInfo.defaultValues = [1 100 100 ones(1,length(networkInfo.sweepVarList)-3)];
    networkInfo.sweepFlags = [0 0 0 ones(1,length(networkInfo.sweepVarList)-3)]==1;                      
    networkInfo.fourStateInputFlags = [1 0 0 ones(1,length(networkInfo.sweepVarList)-3)]==1;
    networkInfo.fourStateWrongInputFlags = [0 1 1 ones(1,length(networkInfo.sweepVarList)-3)]==1;
    networkInfo.paramBounds = repmat([-4 ; 4 ], 1, length(networkInfo.sweepFlags));    
    
    % apply specific constraints
    gammaaFlags = strcmp(networkInfo.sweepVarStrings,'g');    
    networkInfo.paramBounds(1,gammaaFlags) = 0; % gammaa>=1
    deactivationFlags = contains(networkInfo.sweepVarStrings,{'wa','wi'});
    networkInfo.paramBounds(1,deactivationFlags) = 0;    
    
    % networkInfo.paramBounds(1,end-3:end) = 0;
    networkInfo.cr_index = 1;
    networkInfo.cw_index = 2;
    networkInfo.b_index = 3;

    % Provide info for equilibrium constraint application and cycle flux
    % calculations
    networkInfo.fluxFactorFlags = contains(networkInfo.sweepVarStrings,'g');        

    % HM constraints
    networkInfo.hmRatePairs = [ka; ki; wa]; 
    networkInfo.hmRatePairIndices = NaN(size(networkInfo.hmRatePairs));
    for i = 1:size(networkInfo.hmRatePairIndices,1)
        networkInfo.hmRatePairIndices(i,:) = find(ismember(networkInfo.sweepVarList,networkInfo.hmRatePairs(i,:)));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform basic assignments based on type
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    RSymRaw = sym(zeros(nStates));
    wp = 1;
%     wa = 1;
    wi = wa;
    % basic binding and unbinding rates (right and wrong factors have same base
    % rate)
    RSymRaw(ismember(transitionInfoArray,[1,2]) & ~activity_vec_full) = kp/wp;
    RSymRaw(ismember(transitionInfoArray,[1,2]) & activity_vec_full) = kp;

    RSymRaw(ismember(transitionInfoArray,-[1,2])) = km;
%     RSym(ismember(transitionInfoArray,-[1,2]) & activity_vec_full) = kma;

    % basic locus fluctuation rates    
    RSymRaw(ismember(transitionInfoArray,3) & n_total_bound==0) = ka;
    RSymRaw(ismember(transitionInfoArray,3) & n_total_bound==1) = ka*wa*g;    

    RSymRaw(ismember(transitionInfoArray,-3) & n_total_bound==0) = ki;
    RSymRaw(ismember(transitionInfoArray,-3) & n_total_bound==1) = ki/wi;    

    % add specificity and concentration factors
    RSymRaw(transitionInfoArray==1) = RSymRaw(transitionInfoArray==1)*cr;
    RSymRaw(transitionInfoArray==2) = RSymRaw(transitionInfoArray==2)*cw;
    RSymRaw(transitionInfoArray==-2) = RSymRaw(transitionInfoArray==-2)*a;   

    % add diagonal factors 
    RSym = RSymRaw;
    RSym(eye(size(RSymRaw,1))==1) = -sum(RSymRaw);
    
    % construct right and wrong graphs from main fraph
    RSymRight = RSymRaw(1:4,1:4);
    RSymRight(eye(size(RSymRight,1))==1) = -sum(RSymRight);
    
    RSymWrong = RSymRaw([1 6 5 4],[1 6 5 4]);
    RSymWrong(eye(size(RSymWrong,1))==1) = -sum(RSymWrong);
    
    % calculate basic occupancy and sharpness metrics
    
    [V,D] = eig(RSymRight);
    DLog = logical(D==0);
    ssInd = find(all(DLog));
    ssVecSymRight = V(:,ssInd) / sum(V(:,ssInd));    

    % same for wrong network
    [VW,DW] = eig(RSymWrong);
    DWLog = logical(DW==0);
    ssIndW = find(all(DWLog));
    ssVecSymWrong = VW(:,ssIndW) / sum(VW(:,ssIndW));    
    % toc

    % full steady state vector 
    matlabFunction(ssVecSymRight,'File',[savePath 'ssVecFunction4StateRight'],'Optimize',true,'Vars',...
             networkInfo.sweepVarList(networkInfo.fourStateInputFlags));

    % production rate
    productionRateSymRight = sum(ssVecSymRight(activity_vec_full(1:4)));
    matlabFunction(productionRateSymRight,'File',[savePath 'productionRateFunctionFourState'],...
            'Optimize',true,'Vars',networkInfo.sweepVarList(networkInfo.fourStateInputFlags));

    % production rate for hypothetical "wrong" network
    productionRateWrongSym = sum(ssVecSymWrong(activity_vec_full(1:4)));
    matlabFunction(productionRateWrongSym,'File',[savePath 'productionRateWrongFunctionFourState'],...
            'Optimize',true,'Vars',networkInfo.sweepVarList(networkInfo.fourStateWrongInputFlags));

%     % take derivative wrpt c to get sharpness
%     sharpnessSym = diff(productionRateSymRight,cr);
%     matlabFunction(sharpnessSym,'File',[savePath 'sharpnessFunction'],'Optimize',true,'Vars',networkInfo.sweepVarList);
    
    % save helper variables
    activeStates = find(activity_vec_full);
    networkInfo.nStates = nStates;
    networkInfo.RSym = RSym;
    networkInfo.RSymRight = RSymRight;
    networkInfo.RSymWrong = RSymWrong;
    networkInfo.activeStates = activeStates;    
    networkInfo.transitionInfoArray = transitionInfoArray;
    networkInfo.n_right_bound = n_right_bound;
    networkInfo.n_wrong_bound = n_wrong_bound;
    networkInfo.n_total_bound = n_total_bound;
    networkInfo.activeStateFilter = activity_vec_full;