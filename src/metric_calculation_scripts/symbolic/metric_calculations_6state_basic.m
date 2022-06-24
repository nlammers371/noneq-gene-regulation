% This script seeks to builds symbolic functions to calculate gene circuit
% performance characteristics for the 6 state system shown in Main Text
% Figure 4B. Note that it generally takes 20-60 minutes to run due to the
% complexity of the symbolic manipulations performed
clear
close all

% generate path to save metric functions 
subfolderName = 'symbolic';
rootPath = handlePathOptions(subfolderName);

% first, generate baseline 6 state network that serves as kernel for all
% architectures
networkInfo = generateBaselineSixStateNetwork;

% generate prefix for function directory 
folder_prefix = 'n006_s01_ns01_g01';
writePath = [rootPath folder_prefix filesep];
mkdir(writePath)

% helper vector
RSym = networkInfo.RSym;

% designate parameters for 4 state network
nStates = 6;
stateIndex = 1:nStates;
activeStates = [3 4 5];

% RSym = RSym;
% RSym = RSym6(1:nStates,1:nStates);
% % RSymWrong = RSymFull([1 6 5 4],[1 6 5 4]);
% 
% % add diagonal factors 
% RSym(eye(size(RSym,1))==1) = 0;
% RSym(eye(size(RSym,1))==1) = -sum(RSym);
% RSymWrong(eye(size(RSymWrong,1))==1) = -sum(RSymWrong);

% save
RSymFun = matlabFunction(RSym,'File',[writePath 'RSymFun'],'Optimize',true,...
          'Vars',networkInfo.sweepVarList);   
        
% RSymFunWrong = matlabFunction(RSymWrong,'File',[writePath 'RSymFunWrong'],'Optimize',true,...
%           'Vars',networkInfo.sweepVarListWrong);           

% save helper variables
networkInfo.nStates = nStates;
networkInfo.activeStates = activeStates;
networkInfo.permittedConnections = networkInfo.permittedConnections(1:nStates,1:nStates);
networkInfo.transitionInfoArray = networkInfo.transitionInfoArray(1:nStates,1:nStates);
networkInfo.n_right_bound = networkInfo.n_right_bound(1:nStates);
networkInfo.n_wrong_bound = networkInfo.n_wrong_bound(1:nStates);
networkInfo.n_total_bound = networkInfo.n_total_bound(1:nStates);
networkInfo.activeStateFilter = networkInfo.activeStateFilter(1:nStates);

% update key fields
keepVarFlags = true(1,11);
networkInfo.sweepVarList = networkInfo.sweepVarList(keepVarFlags);
networkInfo.sweepVarStrings = networkInfo.sweepVarStrings(keepVarFlags);
networkInfo.defaultValues = networkInfo.defaultValues(keepVarFlags);
networkInfo.sweepFlags = networkInfo.sweepFlags(keepVarFlags);
networkInfo.paramBounds = networkInfo.paramBounds(keepVarFlags);
% networkInfo = rmfield(networkInfo,{'a_index','cw_index'});
networkInfo.forwardRateIndices{1} = networkInfo.forwardRateIndices{1};
networkInfo.backwardRateIndices{1} = networkInfo.backwardRateIndices{1};
networkInfo.RSym = RSym;
save([writePath 'networkInfo'],'networkInfo');       
        
%% %%%%%%%%% derive expressions for production rate and sharpness %%%%%%%%%%

ssVecSym = sixStateSSFromMathematica_v2;

% [V,D] = eig(RSym);
% DLog = logical(D==0);
% ssInd = find(all(DLog));
% ssVecSym = V(:,ssInd) / sum(V(:,ssInd));
% ssVecSym = ssVecSym';

% same for wrong network
% [VW,DW] = eig(RSymWrong);
% DWLog = logical(DW==0);
% ssIndW = find(all(DWLog));
% ssVecSymWrong = VW(:,ssIndW) / sum(VW(:,ssIndW));
% ssVecSymWrong = ssVecSymWrong';
% toc
% %%
% for i = 1:length(networkInfo.sweepVarStrings)
%     syms networkInfo.sweepVarStrings{i},'positive');% cr kmi kma kpi kpa kap kam kip kim positive
% end

% initialize core rate multiplier variables
syms cr cw a positive % cw a

% initialize baseline reaction rates
syms ki ka ku kb positive

% initialize interaction terms
syms wab wib wba wua positive

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% full steady state vector 
matlabFunction(ssVecSym,'File',[writePath 'steadyStateVecFunction'],'Optimize',true,'Vars',networkInfo.sweepVarList);

% production rate
productionRateSym = sum(ssVecSym(activeStates));
productionRateFun = matlabFunction(productionRateSym,'File',[writePath 'productionRateFunction'],'Optimize',true,'Vars',networkInfo.sweepVarList);

% production rate for hypothetical "wrong" network
% productionRateWrongSym = sum(ssVecSymWrong(activeStates));
% productionRateWrongFun = matlabFunction(productionRateWrongSym,'File',[writePath 'productionRateWrongFunction'],'Optimize',true,'Vars',networkInfo.sweepVarListWrong);

% take derivative wrpt c to get sharpness
sharpnessSym = diff(productionRateSym,cr);
sharpnessFun = matlabFunction(sharpnessSym,'File',[writePath 'sharpnessFunction'],'Optimize',true,'Vars',networkInfo.sweepVarList);

% take derivative wrpt c for correct state only
% sharpnessRightSym = diff(ssVecSym(activeStatesFull(1:2)),cr);
% matlabFunction(sharpnessRightSym,'File',[writePath 'sharpnessRightFunction'],'Optimize',true,'Vars',networkInfo.sweepVarList);

%%% %%%%%%%%%%%%%%%%%%%%%%%% Half max constraints %%%%%%%%%%%%%%%%%%%%%%%%%%
% for the 4 state case we can do everything in matlab 
% hmPath = [writePath 'HM_Functions' filesep];
% mkdir(hmPath);
% 
% rate_pairs = [kpi, kpa; kmi, kma; kim, kip; kam kap];
% % rate_index_pairs = [3, 6; 1, 8; 2, 4; 7, 5];
% hmStruct = struct;
% syms a positive;
% 
% for r = 1:size(rate_pairs,1)
%     tic
%     % construct equation to be solved
%     hmSys = productionRateSym == 0.5;
%     hmSys = subs(hmSys,rate_pairs(r,1),rate_pairs(r,1)*a);
%     hmSys = subs(hmSys,rate_pairs(r,2),rate_pairs(r,2)*a);
%     hmSol = solve(hmSys,a,'ReturnConditions',true);    
%     
%     fNames = fieldnames(hmSol);
%     for j = 1:2
%         matlabFunction(hmSol.a(j),'File',[hmPath 'hmFun' num2str(r) num2str(j)],'Optimize',true,'Vars',networkInfo.sweepVarList);
%         matlabFunction(hmSol.conditions(j),'File',[hmPath 'conditionsFun' num2str(r) num2str(j)],'Optimize',true,'Vars',networkInfo.sweepVarList);
%     end
%     hmStruct(r).output_pair = rate_pairs(r,:);
%     hmStruct(r).output_indices = find(ismember(networkInfo.sweepVarList,rate_pairs(r,:)));
%     toc
% end
% 
% save([hmPath 'hmStruct.mat'],'hmStruct');

%%
% generate and solve systems of equations for expected first passage time
% to each active state
activeStateFilter = networkInfo.activeStateFilter;

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
f_vec(activeStateFilter) = 1; % assume all active states produce at the same rate

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
VarFun = matlabFunction(varSum,'File',[writePath 'intrinsicVarianceFunction'],'Optimize',true,'Vars',networkInfo.sweepVarList);
toc

%% %%%%%%%%%%%%%%%%%%%%%%%%% Cycle time calculations %%%%%%%%%%%%%%%%%%%%%%%

%%% mean time to go OFF->ON %%%

% create symbolic vector of passage times
ETVecON = sym('ET%dON', [1 nStates]);
ETVecON(activeStateFilter) = 0;

% create adjusted transition matrix
RSymON = RSym;
Rdiag = -reshape(diag(RSymON),1,[]);
RSymON = RSymON ./ Rdiag;
RSymON(:,activeStateFilter) = 0;
RSymON(eye(size(RSymON))==1) = 0;

% generate system of equations and solve
eqSysON = ETVecON * RSymON + 1./Rdiag;
offStateFilter = ~activeStateFilter;
eqSysON = eqSysON(offStateFilter);
eqSysON = eqSysON == ETVecON(offStateFilter);
eqSolON = solve(eqSysON,ETVecON(offStateFilter));

% transform results into vector and calculate the weighted avers
solVecON = struct2array(eqSolON);

% calculate inbound flux into each OFF state from the ON States
inFluxVecOFF = RSym(offStateFilter,~offStateFilter) * ssVecSym(~offStateFilter)';
ETONMean = (solVecON*inFluxVecOFF) / sum(inFluxVecOFF);
ETONFun = matlabFunction(ETONMean,'File',[writePath 'TauOFFFunction'],'Optimize',true,'Vars',networkInfo.sweepVarList);

%%% mean time to go ON->OFF %%%
onStateFilter = activeStateFilter;

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
ETOFFFun = matlabFunction(ETOFFMean,'File',[writePath 'TauONFunction'],'Optimize',true,'Vars',networkInfo.sweepVarList);

%%% Cycle Time %%%
TauCycleFun = matlabFunction(ETONMean+ETOFFMean,'File',[writePath 'TauCycleFunction'],'Optimize',true,'Vars',networkInfo.sweepVarList);

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

matlabFunction(entropyRateSym,'File',[writePath 'entropyRateFunction'],'Optimize',true,'Vars',networkInfo.sweepVarList);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Copy functions to subdirectories in other folders
% functionList = dir([writePath '*.m']);
% suffix = 'FourState';
% 
% for c = 1:length(copyPaths)
%     mkdir(copyPaths{c});
%     for f = 1:length(functionList)
%         rawName = functionList(f).name;
%         newName = [rawName(1:end-2) suffix '.m'];
%         copyfile([writePath rawName], [copyPaths{c} newName]);
%     end
% end    