% This script seeks to recapitulate calculations previously undertaken in
% mathematica natively in matlab
clear
close all

% generate path to save metric functions 
savePath = '../utilities/metricFunctions/n4_OR/';
copyPaths = {'../utilities/metricFunctions/n6_OR/n4_OR/','../utilities/metricFunctions/n6_OR_NUM/n4_OR/'...
              '../utilities/metricFunctions/n18_OR_NUM/n4_OR/'};
mkdir(savePath);

addpath(genpath('../'))
rmpath(genpath('../utilities/metricFunctions/'));
addpath(savePath);


% define some basic parameters
activeStatesFull = [3 4 5];
activeStates = [3 4];
baseNum = 6;
nStates = 4;

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
activity_vec_full(activeStatesFull) = 1;

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
networkInfo.sweepVarList = [cr kip kap kim kam kpi kmi kpa kma];
networkInfo.sweepVarListWrong = [cw b kip kap kim kam kpi kmi kpa kma];
networkInfo.defaultValues = [1 1 1 1 1 1 1 1 1];
networkInfo.sweepFlags = [0  1 1 1 1 1 1 1 1]==1;
networkInfo.bindingFlags = [0 0 0 0 0 1 0 1 0]==1;
networkInfo.unbindingFlags = [0 0 0 0 0 0 1 0 1]==1;
networkInfo.paramBounds = repmat([-4 ; 4 ], 1, length(networkInfo.sweepFlags));
networkInfo.paramBounds(:,1) = 0;
networkInfo.cr_index = 1;

% Provide info for equilibrium constraint application and cycle flux
% calculations
networkInfo.forwardRateConstants = [kpi kap kma kim];
networkInfo.forwardRateIndices = {find(ismember(networkInfo.sweepVarList,networkInfo.forwardRateConstants))};
networkInfo.backwardRateConstants = [kmi kip kpa kam];
networkInfo.backwardRateIndices = {find(ismember(networkInfo.sweepVarList,networkInfo.backwardRateConstants))};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform basic assignments based on type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RSym = sym(zeros(baseNum));

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

RSymFull = RSym;
RSym = RSym(1:nStates,1:nStates);
RSymWrong = RSymFull([1 6 5 4],[1 6 5 4]);

% add diagonal factors 
RSym(eye(size(RSym,1))==1) = -sum(RSym);
RSymWrong(eye(size(RSymWrong,1))==1) = -sum(RSymWrong);

% save
RSymFun = matlabFunction(RSym,'File',[savePath 'RSymFun'],'Optimize',true,...
          'Vars',networkInfo.sweepVarList);   
        
RSymFunWrong = matlabFunction(RSymWrong,'File',[savePath 'RSymFunWrong'],'Optimize',true,...
          'Vars',networkInfo.sweepVarListWrong);           

% save helper variables
networkInfo.nStates = nStates;
networkInfo.activeStates = activeStatesFull;
networkInfo.permittedConnections = permittedConnections(1:nStates,1:nStates);
networkInfo.transitionInfoArray = transitionInfoArray(1:nStates,1:nStates);
networkInfo.n_right_bound = n_right_bound(1:nStates);
networkInfo.n_wrong_bound = n_wrong_bound(1:nStates);
networkInfo.n_total_bound = n_total_bound(1:nStates);
networkInfo.activeStateFilter = activity_vec_full(1:nStates);

save([savePath 'networkInfo'],'networkInfo');       
        
%% %%%%%%%%% derive expressions for production rate and sharpness %%%%%%%%%%

[V,D] = eig(RSym);
DLog = logical(D==0);
ssInd = find(all(DLog));
ssVecSym = V(:,ssInd) / sum(V(:,ssInd));
ssVecSym = ssVecSym';

% same for wrong network
[VW,DW] = eig(RSymWrong);
DWLog = logical(DW==0);
ssIndW = find(all(DWLog));
ssVecSymWrong = VW(:,ssIndW) / sum(VW(:,ssIndW));
ssVecSymWrong = ssVecSymWrong';
% toc
% %%
syms cr kmi kma kpi kpa kap kam kip kim positive

% full steady state vector 
matlabFunction(ssVecSym,'File',[savePath 'steadyStateVecFunction'],'Optimize',true,'Vars',networkInfo.sweepVarList);

% production rate
productionRateSym = sum(ssVecSym(activeStates));
productionRateFun = matlabFunction(productionRateSym,'File',[savePath 'productionRateFunction'],'Optimize',true,'Vars',networkInfo.sweepVarList);

% production rate for hypothetical "wrong" network
productionRateWrongSym = sum(ssVecSymWrong(activeStates));
productionRateWrongFun = matlabFunction(productionRateWrongSym,'File',[savePath 'productionRateWrongFunction'],'Optimize',true,'Vars',networkInfo.sweepVarListWrong);

% take derivative wrpt c to get sharpness
sharpnessSym = diff(productionRateSym,cr);
sharpnessFun = matlabFunction(sharpnessSym,'File',[savePath 'sharpnessFunction'],'Optimize',true,'Vars',networkInfo.sweepVarList);

% take derivative wrpt c for correct state only
sharpnessRightSym = diff(ssVecSym(activeStatesFull(1:2)),cr);
matlabFunction(sharpnessRightSym,'File',[savePath 'sharpnessRightFunction'],'Optimize',true,'Vars',networkInfo.sweepVarList);

%% %%%%%%%%%%%%%%%%%%%%%%%% Half max constraints %%%%%%%%%%%%%%%%%%%%%%%%%%
% for the 4 state case we can do everything in matlab 
hmPath = [savePath 'HM_Functions' filesep];
mkdir(hmPath);

rate_pairs = [kpi, kpa; kmi, kma; kim, kip; kam kap];
% rate_index_pairs = [3, 6; 1, 8; 2, 4; 7, 5];
hmStruct = struct;
syms a positive;

for r = 1:size(rate_pairs,1)
    tic
    % construct equation to be solved
    hmSys = productionRateSym == 0.5;
    hmSys = subs(hmSys,rate_pairs(r,1),rate_pairs(r,1)*a);
    hmSys = subs(hmSys,rate_pairs(r,2),rate_pairs(r,2)*a);
    hmSol = solve(hmSys,a,'ReturnConditions',true);    
    
    fNames = fieldnames(hmSol);
    for j = 1:2
        matlabFunction(hmSol.a(j),'File',[hmPath 'hmFun' num2str(r) num2str(j)],'Optimize',true,'Vars',networkInfo.sweepVarList);
        matlabFunction(hmSol.conditions(j),'File',[hmPath 'conditionsFun' num2str(r) num2str(j)],'Optimize',true,'Vars',networkInfo.sweepVarList);
    end
    hmStruct(r).output_pair = rate_pairs(r,:);
    hmStruct(r).output_indices = find(ismember(networkInfo.sweepVarList,rate_pairs(r,:)));
    toc
end

save([hmPath 'hmStruct.mat'],'hmStruct');

%%
% generate and solve systems of equations for expected first passage time
% to each active state

etInfo = struct;
for CurrentState = 1:nStates
      
    % create symbolic vector of passage times
    etInfo(CurrentState).ETVec = sym(['ET%d' num2str(CurrentState)], [1 nStates]);
    etInfo(CurrentState).ETVec(CurrentState) = 0;
    
    % create adjusted transition matrix
    etInfo(CurrentState).RSymTo = RSym;
    etInfo(CurrentState).Rdiag = -reshape(diag(etInfo(CurrentState).RSymTo),1,[]);
    etInfo(CurrentState).RSymTo = etInfo(CurrentState).RSymTo ./ etInfo(CurrentState).Rdiag;
    etInfo(CurrentState).RSymTo(:,CurrentState) = 0;
    etInfo(CurrentState).RSymTo(eye(size(etInfo(CurrentState).RSymTo))==1) = 0;
    
    % create system of equations
    etInfo(CurrentState).eqSys = etInfo(CurrentState).ETVec * etInfo(CurrentState).RSymTo + 1./etInfo(CurrentState).Rdiag;
    
    stateFilter = stateIndex~=CurrentState;    
    
    etInfo(CurrentState).eqSys = etInfo(CurrentState).eqSys(stateFilter);
    etInfo(CurrentState).eqSys = etInfo(CurrentState).eqSys == etInfo(CurrentState).ETVec(stateFilter);
    
    % solve
    etInfo(CurrentState).eqSol = solve(etInfo(CurrentState).eqSys,etInfo(CurrentState).ETVec(stateFilter));
    etInfo(CurrentState).etSolVec = struct2array(etInfo(CurrentState).eqSol);  
    etInfo(CurrentState).etSolVec = [etInfo(CurrentState).etSolVec(1:CurrentState-1) 0 etInfo(CurrentState).etSolVec(CurrentState:end)];
    etInfo(CurrentState).ETMean = etInfo(CurrentState).etSolVec*ssVecSym';
end
 
% %%%%%%%%%%%%%%%%%%% get expressions for variance %%%%%%%%%%%%%%%%%%%%%%%


f_vec = zeros(1,nStates); % initiation rate for each state
f_vec(activeStatesFull) = 1; % assume all active states produce at the same rate

% construct Z matrix 
% see eqs 28 and 29 from: "Asymptotic Formulas for Markov Processes with
%                          Applications to Simulation"
ZSym = sym(size(RSym));
for i = 1:nStates
    for j = 1:nStates
        if i == j
            ZSym(i,j) = ssVecSym(j)*etInfo(j).ETMean;
        else
            ZSym(i,j) = ssVecSym(j)*(etInfo(j).ETMean-etInfo(j).etSolVec(i));
        end
    end
end

% now caculate the variance (see eq 12 from above source)
varSum = 0;
for i = 1:size(RSym,2)
    for j = 1:size(RSym,2)
        varSum = varSum + f_vec(i)*ssVecSym(i)*ZSym(i,j)*f_vec(j);
    end
end
varSum = 2*varSum;

% convert to function
tic
VarFun = matlabFunction(varSum,'File',[savePath 'intrinsicVarianceFunction'],'Optimize',true,'Vars',networkInfo.sweepVarList);
toc
%% %%%%%%%%%%%%%%%%%%%%%%%%% Cycle time calculations %%%%%%%%%%%%%%%%%%%%%%%

%%% mean time to go OFF->ON %%%

% create symbolic vector of passage times
ETVecON = sym('ET%dON', [1 nStates]);
ETVecON(activeStatesFull) = 0;

% create adjusted transition matrix
RSymON = RSym;
Rdiag = -reshape(diag(RSymON),1,[]);
RSymON = RSymON ./ Rdiag;
RSymON(:,activeStatesFull) = 0;
RSymON(eye(size(RSymON))==1) = 0;

% generate system of equations and solve
eqSysON = ETVecON * RSymON + 1./Rdiag;
offStateFilter = ~ismember(stateIndex,activeStatesFull);
eqSysON = eqSysON(offStateFilter);
eqSysON = eqSysON == ETVecON(offStateFilter);
eqSolON = solve(eqSysON,ETVecON(offStateFilter));

% transform results into vector and calculate the weighted avers
solVecON = struct2array(eqSolON);

% calculate inbound flux into each OFF state from the ON States
inFluxVecOFF = RSym(offStateFilter,~offStateFilter) * ssVecSym(~offStateFilter)';
ETONMean = (solVecON*inFluxVecOFF) / sum(inFluxVecOFF);
ETONFun = matlabFunction(ETONMean,'File',[savePath 'TauOFFFunction'],'Optimize',true,'Vars',networkInfo.sweepVarList);

%%% mean time to go ON->OFF %%%
onStateFilter = ismember(stateIndex,activeStatesFull);

% create symbolic vector of passage times
ETVecOFF = sym('ET%dOFF', [1 nStates]);
ETVecOFF(~onStateFilter) = 0;

% create adjusted transition matrix
RSymOFF = RSym;
Rdiag = -reshape(diag(RSymOFF),1,[]);
RSymOFF = RSymOFF ./ Rdiag;

RSymOFF(:,~onStateFilter) = 0;
RSymOFF(eye(size(RSymOFF))==1) = 0;

% generate system of equations and solve
eqSysOFF = ETVecOFF * RSymOFF + 1./Rdiag;
eqSysOFF = eqSysOFF(onStateFilter);
eqSysOFF = eqSysOFF == ETVecOFF(onStateFilter);
eqSolOFF = solve(eqSysOFF,ETVecOFF(onStateFilter));

% transform results into vector and calculate the weighted avers
solVecOFF = struct2array(eqSolOFF);
% calculate incoming flux from OFF to ON states
inFluxVecON = RSym(onStateFilter,~onStateFilter) * ssVecSym(~onStateFilter)';
ETOFFMean = (solVecOFF*inFluxVecON) / sum(inFluxVecON);
ETOFFFun = matlabFunction(ETOFFMean,'File',[savePath 'TauONFunction'],'Optimize',true,'Vars',networkInfo.sweepVarList);

%%% Cycle Time %%%
TauCycleFun = matlabFunction(ETONMean+ETOFFMean,'File',[savePath 'TauCycleFunction'],'Optimize',true,'Vars',networkInfo.sweepVarList);

%%% %%%%%%%%%%%%%%%%%%%% Energy Dissipation Expressions %%%%%%%%%%%%%%%%%%%%

%%% Entropy production %%%
% See equation 5 in: Thermodynamics of Statistical Inference by Cells

entropyRateSym = 0;

for i = 1:nStates
    for j = 1:nStates
        if i ~= j && RSym(i,j)~=0
            entropyRateSym = entropyRateSym + ssVecSym(i)*RSym(j,i)*log(RSym(j,i)/RSym(i,j));
        end
    end
end

matlabFunction(entropyRateSym,'File',[savePath 'entropyRateFunction'],'Optimize',true,'Vars',networkInfo.sweepVarList);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Copy functions to subdirectories in other folders
functionList = dir([savePath '*.m']);
suffix = 'FourState';

for c = 1:length(copyPaths)
    mkdir(copyPaths{c});
    for f = 1:length(functionList)
        rawName = functionList(f).name;
        newName = [rawName(1:end-2) suffix '.m'];
        copyfile([savePath rawName], [copyPaths{c} newName]);
    end
end    