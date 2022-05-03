% This script explores possibility of keeping SS vec in symbolic form
clear
close all

savePath = '../utilities/metricFunctions/n12_OR_enhancer_NUM/';
mkdir(savePath);
addpath(savePath)

% define some basic parameters
activeStatesBase = [3 4 5];
baseNum = 6;
nStates = 12;
% sensingEdges = [2,1 ; 3,4];

% define symbolic variables
RSymFull = sym('k%d%d', [nStates nStates],'positive');

% zero out forbidden transitions
[fromRef,toRef] = meshgrid(1:baseNum,1:baseNum);
diffRef = abs(fromRef-toRef);
toRefHalved = toRef<=baseNum/2;
fromRefHalved = fromRef<=baseNum/2;
permittedConnectionsRaw= (diffRef==1 & toRefHalved==fromRefHalved) | diffRef==baseNum/2;

% permute these connections to follow a more intuitive labeling scheme
indexCurr = 1:baseNum;
indexAdjusted = circshift(indexCurr,floor(baseNum/4));
indexAdjusted = [indexAdjusted(1:baseNum/2) fliplr(indexAdjusted(baseNum/2+1:end))];
[~,si] = sort(indexAdjusted);
permittedConnectionsInit = permittedConnectionsRaw(si,si);

% generate an array with binding info
transitionInfoInit = zeros(size(permittedConnectionsInit));

% specific binding/unbinding
spec_pairs = {[1,2],[4,3]};
for s = 1:length(spec_pairs)
    ind1 = spec_pairs{s}(1);
    ind2 = spec_pairs{s}(2);
    % update
    transitionInfoInit(ind2,ind1) = 1;
    transitionInfoInit(ind1,ind2) = -1;
end

% non-specific binding/unbinding
non_spec_pairs = {[1,6],[4,5]};
for s = 1:length(non_spec_pairs)
    ind1 = non_spec_pairs{s}(1);
    ind2 = non_spec_pairs{s}(2);
    % update
    transitionInfoInit(ind2,ind1) = 2;
    transitionInfoInit(ind1,ind2) = -2;
end

% locus activity fluctuations
locus_pairs = {[1,4],[6,5],[2,3]};
for s = 1:length(locus_pairs)
    ind1 = locus_pairs{s}(1);
    ind2 = locus_pairs{s}(2);
    % update
    transitionInfoInit(ind2,ind1) = 3;
    transitionInfoInit(ind1,ind2) = -3;
end

% generate array that indicates activity state
activity_vec_init = false(1,baseNum);
activity_vec_init(activeStatesBase) = 1;

% generate flags indicating numbers of right and wrong factors bound
n_right_bound_init = zeros(size(activity_vec_init));
n_wrong_bound_init = n_right_bound_init;
n_right_bound_init([2 3]) = 1;
n_wrong_bound_init([5 6]) = 1;

% generate full 12 state array
rate_mask = repelem(eye(2),baseNum,baseNum);

% full connection array
permittedConnections = repmat(permittedConnectionsInit,2,2);%permittedConnectionsInit
permittedConnections = permittedConnections.*rate_mask;

% full binding array
transitionInfoArray = repmat(transitionInfoInit,2,2);%permittedConnectionsInit
transitionInfoArray = transitionInfoArray.*rate_mask;

% full activity vector
activity_vec_full = repmat(activity_vec_init,1,2);

% binding info vecs
n_right_bound = repmat(n_right_bound_init,1,2);

n_wrong_bound = repmat(n_wrong_bound_init,1,2);

n_total_bound = n_wrong_bound + n_right_bound;

% create mask indicating enhancer association state
enhancer_associated_flags = false(size(transitionInfoArray));
enhancer_associated_flags(:,baseNum+1:end) = true;

% add cross-plane connections. Let us assume the 6 state plane that is
% equivalent to the simpler 6 state model considered correspond to the
% first block
for i = 1:baseNum    
    for j = 1
        ind1 = i;
        ind2 = baseNum*j + i;
        
        % add to array
        permittedConnections(ind1,ind2) = 1;
        permittedConnections(ind2,ind1) = 1;
        
        % add binding info
        transitionInfoArray(ind2,ind1) = 4; % specify whether it is a cognate or non-cognate factor binding
        transitionInfoArray(ind1,ind2) = -4; % all unbinding events from non-specific site are equivalent
    end
end
permittedConnections = permittedConnections==1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate symbolic transition rate matrix

% initialize core rate multiplier variables
syms cr cw b positive

% initialize core locus activity transition rates
syms kip kap kim kam positive

% initialize core binding rates
syms kpi kpa kmi kma positive

% initialize enhancer-related rates
syms kem kep kdp kdm positive % assume 1 and only 1 association rate, but two different disassociation rates depending on binding status

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specify var order, var bounds, and whether it can be swept by default
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
networkInfo = struct;
networkInfo.sweepVarList = [cr cw b kip kap kim kam kpi kmi kpa kma kep kem kdp kdm];
networkInfo.defaultValues = [1 100 100 1 1 1 1 1 1 1 1 1 1 1 1];
networkInfo.sweepFlags = [0 0 0 1 1 1 1 1 1 1 1 1 1 1 1]==1;
networkInfo.fourStateInputFlags = [1 0 0 1 1 1 1 1 1 1 1 0 0 0 0]==1;
networkInfo.fourStateWrongInputFlags = [0 1 1 1 1 1 1 1 1 1 1 0 0 0 0]==1;
networkInfo.bindingFlags = [0 0 0 0 0 0 0 1 0 1 0 0 0 0 0];
networkInfo.unbindingFlags = [0 0 0 0 0 0 0 0 1 0 1 0 0 0 0];
networkInfo.paramBounds = repmat([-4 ; 4 ], 1, length(networkInfo.sweepFlags));
networkInfo.paramBounds(:,1:3) = 0;
% networkInfo.paramBounds(1,end-3:end) = 0;
networkInfo.cr_index = 1;
networkInfo.cw_index = 2;
networkInfo.b_index = 3;

% Provide info for equilibrium constraint application and cycle flux
% calculations

% for 12 state network, we have to start with cooperativity and binding
% factors
networkInfo.forwardRateConstants(1) = {[kpi kma]};
networkInfo.forwardRateIndices(1) = {find(ismember(networkInfo.sweepVarList,networkInfo.forwardRateConstants{1}))};
networkInfo.backwardRateConstants(1) = {[kpa kmi]};
networkInfo.backwardRateIndices(1) = {find(ismember(networkInfo.sweepVarList,networkInfo.backwardRateConstants{1}))};
  
networkInfo.forwardRateConstants(2) = {[kap kim]};
networkInfo.forwardRateIndices(2) = {find(ismember(networkInfo.sweepVarList,networkInfo.forwardRateConstants{2}))};
networkInfo.backwardRateConstants(2) = {[ kip kam]};
networkInfo.backwardRateIndices(2) = {find(ismember(networkInfo.sweepVarList,networkInfo.backwardRateConstants{2}))};

networkInfo.forwardRateConstants(3) = {[kdm kep]};
networkInfo.forwardRateIndices(3) = {find(ismember(networkInfo.sweepVarList,networkInfo.forwardRateConstants{3}))};
networkInfo.backwardRateConstants(3) = {[kdp kem]};
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
% create some useful filter variables


% basic binding and unbinding rates (right and wrong factors have same base
% rate)
RSym(ismember(transitionInfoArray,[1,2]) & ~activity_vec_full) = kpi;
RSym(ismember(transitionInfoArray,[1,2]) & activity_vec_full) = kpa;

RSym(ismember(transitionInfoArray,-[1,2]) & ~activity_vec_full) = kmi;
RSym(ismember(transitionInfoArray,-[1,2]) & activity_vec_full) = kma;

% basic locus fluctuation rates
% assume that the number bound (but not the identity) can alter this rate
RSym(ismember(transitionInfoArray,3) & n_total_bound==0) = kam;
RSym(ismember(transitionInfoArray,3) & n_total_bound==1 & enhancer_associated_flags) = kap;
RSym(ismember(transitionInfoArray,3) & n_total_bound==1 & ~enhancer_associated_flags) = kam; % it doesn't matter if the activator is bound but enhancer not associated

RSym(ismember(transitionInfoArray,-3) & n_total_bound==0) = kim;
RSym(ismember(transitionInfoArray,-3) & n_total_bound==1 & enhancer_associated_flags) = kip;
RSym(ismember(transitionInfoArray,-3) & n_total_bound==1 & ~enhancer_associated_flags) = kim;

% enhancer dynamics
RSym(ismember(transitionInfoArray,4 & n_total_bound==0)) = kem;
RSym(ismember(transitionInfoArray,4 & n_total_bound==1)) = kep;
% RSym(ismember(transitionInfoArray,3) & n_total_bound==1) = kap;

RSym(ismember(transitionInfoArray,-4) & n_total_bound==0) = kdm;
RSym(ismember(transitionInfoArray,-4) & n_total_bound==1) = kdp;


% add specificity and concentration factors
RSym(transitionInfoArray==1) = RSym(transitionInfoArray==1)*cr;
RSym(transitionInfoArray==2) = RSym(transitionInfoArray==2)*cw;
RSym(transitionInfoArray==-2) = RSym(transitionInfoArray==-2)*b;

% add diagonal factors 
RSym(eye(size(RSym,1))==1) = -sum(RSym);

% wrte to file
RSymFun = matlabFunction(RSym,'File',[savePath 'RSymFun'],'Optimize',true,...
          'Vars',networkInfo.sweepVarList);       
        
% save helper variables
activeStates = find(activity_vec_full);
networkInfo.nStates = nStates;
networkInfo.activeStates = activeStates;
networkInfo.permittedConnections = permittedConnections;
networkInfo.transitionInfoArray = transitionInfoArray;
networkInfo.n_right_bound = n_right_bound;
networkInfo.n_wrong_bound = n_wrong_bound;
networkInfo.n_total_bound = n_total_bound;
networkInfo.activeStateFilter = activity_vec_full;

save([savePath 'networkInfo'],'networkInfo');        
