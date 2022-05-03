% This script generates helper functions for numeric metric calculations
clear
close all

% generate path to save metric functions 
savePath = '../utilities/metricFunctions/n6_OR_NUM/';
mkdir(savePath);
addpath(savePath)

% define some basic parameters
activeStates = [3 4 5];
baseNum = 6;
nStates = 6;

% helper vector
stateIndex = 1:nStates;

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
permittedConnections = permittedConnectionsRaw(si,si);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate an array with binding info
transitionInfoArray = zeros(size(permittedConnections));

% specific binding/unbinding
spec_pairs = {[1,2],[4,3]};
for s = 1:length(spec_pairs)
    ind1 = spec_pairs{s}(1);
    ind2 = spec_pairs{s}(2);
    % update
    transitionInfoArray(ind2,ind1) = 1;
    transitionInfoArray(ind1,ind2) = -1;
end

% non-specific binding/unbinding
non_spec_pairs = {[1,6],[4,5]};
for s = 1:length(non_spec_pairs)
    ind1 = non_spec_pairs{s}(1);
    ind2 = non_spec_pairs{s}(2);
    % update
    transitionInfoArray(ind2,ind1) = 2;
    transitionInfoArray(ind1,ind2) = -2;
end

% locus activity fluctuations
locus_pairs = {[1,4],[6,5],[2,3]};
for s = 1:length(locus_pairs)
    ind1 = locus_pairs{s}(1);
    ind2 = locus_pairs{s}(2);
    % update
    transitionInfoArray(ind2,ind1) = 3;
    transitionInfoArray(ind1,ind2) = -3;
end

% generate array tht indicates activity state
activity_vec_full = false(1,baseNum);
activity_vec_full(activeStates) = 1;

% generate flags indicating numbers of right and wrong factors bound
n_right_bound = zeros(size(activity_vec_full));
n_wrong_bound = n_right_bound;
n_right_bound([2 3]) = 1;
n_wrong_bound([5 6]) = 1;

permittedConnections = permittedConnections==1;
n_total_bound = n_wrong_bound + n_right_bound;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate symbolic transition rate matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize core rate multiplier variables
syms cr cw b positive

% initialize core locus activity transition rates
syms kip kap kim kam positive

% initialize core binding rates
syms kpi kpa kmi kma positive

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specify var order, var bounds, and whether it can be swept by default
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
networkInfo = struct;
networkInfo.sweepVarList = [cr cw b kip kap kim kam kpi kmi kpa kma];
networkInfo.defaultValues = [1 100 100 1 1 1 1 1 1 1 1];
networkInfo.sweepFlags = [0 0 0 1 1 1 1 1 1 1 1]==1;

networkInfo.fourStateInputFlags = [1 0 0 1 1 1 1 1 1 1 1]==1;
networkInfo.fourStateWrongInputFlags = [0 1 1 1 1 1 1 1 1 1 1]==1;

networkInfo.bindingFlags = [0 0 0 0 0 0 0 1 0 1 0];
networkInfo.unbindingFlags = [0 0 0 0 0 0 0 0 1 0 1];

networkInfo.paramBounds = repmat([-4 ; 4 ], 1, length(networkInfo.sweepFlags));
networkInfo.paramBounds(:,1:3) = 0;

networkInfo.cr_index = 1;
networkInfo.cw_index = 2;
networkInfo.b_index = 3;

% Provide info for equilibrium constraint application and cycle flux
% calculations
networkInfo.forwardRateConstants = [kpi kap kma kim];
networkInfo.forwardRateIndices = {find(ismember(networkInfo.sweepVarList,networkInfo.forwardRateConstants))};
networkInfo.backwardRateConstants = [kmi kip kpa kam];
networkInfo.backwardRateIndices = {find(ismember(networkInfo.sweepVarList,networkInfo.backwardRateConstants))};

%%% HM constraints
networkInfo.hmRatePairs = [kpi kpa; kmi kma; kip kim; kap kam]; 
networkInfo.hmRatePairIndices = NaN(size(networkInfo.hmRatePairs));
for i = 1:size(networkInfo.hmRatePairIndices,1)
    networkInfo.hmRatePairIndices(i,:) = find(ismember(networkInfo.sweepVarList,networkInfo.hmRatePairs(i,:)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform basic assignments based on type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RSym = sym(zeros(nStates));

% basic binding and unbinding rates
RSym(ismember(transitionInfoArray,[1,2]) & ~activity_vec_full) = kpi;
RSym(ismember(transitionInfoArray,[1,2]) & activity_vec_full) = kpa;

RSym(ismember(transitionInfoArray,-[1,2]) & ~activity_vec_full) = kmi;
RSym(ismember(transitionInfoArray,-[1,2]) & activity_vec_full) = kma;

% basic locus fluctuation rates
RSym(ismember(transitionInfoArray,3) & n_total_bound==0) = kam;
RSym(ismember(transitionInfoArray,3) & n_total_bound==1) = kap;

RSym(ismember(transitionInfoArray,-3) & n_total_bound==0) = kim;
RSym(ismember(transitionInfoArray,-3) & n_total_bound==1) = kip;

% add specificity and concentration factors
RSym(transitionInfoArray==1) = RSym(transitionInfoArray==1)*cr;
RSym(transitionInfoArray==2) = RSym(transitionInfoArray==2)*cw;
RSym(transitionInfoArray==-2) = RSym(transitionInfoArray==-2)*b;

% add diagonal factors 
RSym(eye(size(RSym,1))==1) = -sum(RSym);

% save
RSymFun = matlabFunction(RSym,'File',[savePath 'RSymFun'],'Optimize',true,...
          'Vars',networkInfo.sweepVarList);   

% save helper variables
networkInfo.nStates = nStates;
networkInfo.activeStates = activeStates;
networkInfo.permittedConnections = permittedConnections;
networkInfo.transitionInfoArray = transitionInfoArray;
networkInfo.n_right_bound = n_right_bound;
networkInfo.n_wrong_bound = n_wrong_bound;
networkInfo.n_total_bound = n_total_bound;
networkInfo.activeStateFilter = activity_vec_full;

save([savePath 'networkInfo'],'networkInfo');