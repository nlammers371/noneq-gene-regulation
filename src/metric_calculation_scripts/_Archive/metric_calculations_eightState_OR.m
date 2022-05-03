% This script seeks to recapitulate calculations previously undertaken in
% mathematica natively in matlab
clear
close all
addpath(genpath('../'))
rmpath(genpath('../utilities/metricFunctions/'));
addpath(['../utilities/metricFunctions/n8_OR/']);

% generate path to save metric functions 
savePath = '../utilities/metricFunctions/n8_OR/';
mkdir(savePath);

% define some basic parameters
activeStates = [3 4 5 8];
nStates = 6; % temporary dummy entry
sensingEdges = [2,1 ; 3,4];
stateIndex = 1:nStates;


%%
% zero out forbidden transitions
[fromRef,toRef] = meshgrid(1:nStates,1:nStates);
diffRef = abs(fromRef-toRef);
toRefHalved = toRef<=nStates/2;
fromRefHalved = fromRef<=nStates/2;
permittedConnections = (diffRef==1 & toRefHalved==fromRefHalved) | diffRef==nStates/2;

% permute these connections to follow a more intuitive labeling scheme
indexCurr = 1:nStates;
indexAdjusted = circshift(indexCurr,floor(nStates/4));
indexAdjusted = [indexAdjusted(1:nStates/2) fliplr(indexAdjusted(nStates/2+1:end))];
[~,si] = sort(indexAdjusted);
permittedConnections = permittedConnections(si,si);

permittedConnections(:,7:8) = 0;
permittedConnections(7:8,:) = 0;
% Apply manual corrections 
%%
% To 8 corrections
permittedConnections(8,[3 5]) = 1;

% From 8 corrections
permittedConnections([3 5],8) = 1;

% To 7
permittedConnections(7,[2 6]) = 1;

% From 7
permittedConnections([2 6],7) = 1;


%%

nStates = 8;
% define symbolic variables
RSymFull = sym('k%d%d', [nStates nStates],'positive');

% update transition matrix to reflect network topology
RSymRaw = RSymFull.*permittedConnections;
  
for i = 1:nStates
    for j = 1:nStates
        if i~=j 
            string = char(RSymRaw(i,j));
            if length(string)>1
                syms(string,'positive')
            end
        end
    end
end
        
%incorporate sensing edges   
syms c positive
sensingIndices = sub2ind(size(RSymRaw),sensingEdges(:,1),sensingEdges(:,2));
RSymSimp = RSymRaw;
RSymSimp(sensingIndices) = c*RSymSimp(sensingIndices);

% add diagonal terms
RSymSimp(eye(size(RSymSimp))==1) = -sum(RSymSimp);        
        
%% apply replacement rules
syms cw c b
RSym = subs(RSymSimp,[k65,k56,k61,k16,k54,k45,k27,k38,k58,k67,k76,k72,k83,k85],...
          [k23,k32,cw*k21,b*k12,cw*k34,b*k43, b*k12,b*k43,k43,k12,c*k21,cw*k21,cw*k34,c*k34]);
  
%% experiment with solving for the SS vector
SSVecFull = sym('ss%d', [1 nStates],'positive');

%%
tic
n = size(RSym,1);
P = [zeros(1,n),1]/[RSym-eye(n),ones(n,1)]
toc
%% %%%%%%%%% derive expressions for production rate and sharpness %%%%%%%%%%
% NL: for now, I copy RSym from mathematica and paste the resulting solution
% for the steady state vector back into  matlab
% tic
% [V,D] = eig(RSym);
% DLog = logical(D==0);
% ssInd = find(all(DLog));
% toc
%%
ssVecSym = eightStateSSFromMathematica(RSymRaw);
% rateSymAlt = sum(ssVecSym(activeStates));
%%
% full steady state vector 
tic
matlabFunction(ssVecSym,'File',[savePath 'steadyStateVecFunction'],'Optimize',true);
toc
%%
% production rate
productionRateSym = sum(ssVecSym(activeStates));
productionRateFun = matlabFunction(productionRateSym,'File',[savePath 'productionRateFunction'],'Optimize',true);

% take derivative wrpt c to get sharpness
sharpnessSym = diff(productionRateSym,c);
sharpnessFun = matlabFunction(sharpnessSym,'File',[savePath 'sharpnessFunction'],'Optimize',true);

% take derivative wrpt c for correct state only
sharpnessRightSym = diff(ssVecSym(activeStates(1:2)),c);
matlabFunction(sharpnessRightSym,'File',[savePath 'sharpnessRightFunction'],'Optimize',true);


%% %%%%%%%%%%%%%%%%%%%%%%%% Half max constraints %%%%%%%%%%%%%%%%%%%%%%%%%%
% to start, we've gotta make a simplified version of the full 6 state model
% that assumes symmetry between right and wrong cycles

for i = 4
    tic
    nFun = 3;
    if i == 1
        nFun = 2;
    end
    for j = 1:nFun
        numStr = [num2str(i) num2str(j)];
        try        
            eval(['a' numStr 'Sym = a' numStr 'FunFromMathematica(RSymRaw);'])
            savePathFull = [savePath 'a' numStr 'SymFun'];
            eval(['matlabFunction(a' numStr 'Sym,"File",savePathFull,"Optimize",true);'])
        catch
            disp(['error in ' numStr])
        end
    end
    toc
end

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
f_vec(activeStates) = 1; % assume all active states produce at the same rate

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
VarFun = matlabFunction(varSum,'File',[savePath 'intrinsicVarianceFunction'],'Optimize',true);
toc
%% %%%%%%%%%%%%%%%%%%%%%%%%% Cycle time calculations %%%%%%%%%%%%%%%%%%%%%%%

%%% mean time to go OFF->ON %%%

% create symbolic vector of passage times
ETVecON = sym('ET%dON', [1 nStates]);
ETVecON(activeStates) = 0;

% create adjusted transition matrix
RSymON = RSym;
Rdiag = -reshape(diag(RSymON),1,[]);
RSymON = RSymON ./ Rdiag;
RSymON(:,activeStates) = 0;
RSymON(eye(size(RSymON))==1) = 0;

% generate system of equations and solve
eqSysON = ETVecON * RSymON + 1./Rdiag;
offStateFilter = ~ismember(stateIndex,activeStates);
eqSysON = eqSysON(offStateFilter);
eqSysON = eqSysON == ETVecON(offStateFilter);
eqSolON = solve(eqSysON,ETVecON(offStateFilter));

% transform results into vector and calculate the weighted avers
solVecON = struct2array(eqSolON);

% calculate inbound flux into each OFF state from the ON States
inFluxVecOFF = RSym(offStateFilter,~offStateFilter) * ssVecSym(~offStateFilter)';
ETONMean = (solVecON*inFluxVecOFF) / sum(inFluxVecOFF);
ETONFun = matlabFunction(ETONMean,'File',[savePath 'TauOFFFunction'],'Optimize',true);

%%% mean time to go ON->OFF %%%
onStateFilter = ismember(stateIndex,activeStates);

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
ETOFFFun = matlabFunction(ETOFFMean,'File',[savePath 'TauONFunction'],'Optimize',true);

%%% Cycle Time %%%
TauCycleFun = matlabFunction(ETONMean+ETOFFMean,'File',[savePath 'TauCycleFunction'],'Optimize',true);

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

matlabFunction(entropyRateSym,'File',[savePath 'entropyRateFunction'],'Optimize',true);

%%% calculate flux (only well-defined for single-loop systems)
% netFluxSym = ssVecSym(1)*RSym(4,1) - ssVecSym(4)*RSym(1,4);
% netFluxFun = matlabFunction(netFluxSym);
% matlabFunction(netFluxSym,'File',[savePath 'netFluxFunction4State']);

% %%% forward and backward fluxes
% [Numerator, Denominator] = numden(productionRateSym);
% forwardFluxSym = RSym(1,2)*RSym(2,3)*RSym(3,4)*RSym(4,1) / Denominator;
% forwardFluxFun = matlabFunction(forwardFluxSym);
% matlabFunction(forwardFluxSym,'File',[savePath 'forwardFluxFunction4State']);
% 
% backwardFluxSym = RSym(2,1)*RSym(3,2)*RSym(4,3)*RSym(1,4) / Denominator;
% backwardFluxFun = matlabFunction(backwardFluxSym);
% matlabFunction(backwardFluxSym,'File',[savePath 'backwardFluxFunction4State']);

