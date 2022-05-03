% This script seeks to recapitulate calculations previously undertaken in
% mathematica natively in matlab
clear
close all

% generate path to save metric functions 
savePath = '../utilities/metricFunctions/fourState/';
mkdir(savePath);

% define some basic parameters
activeStates = [3 4];

% define symbolic variables
reverseVec = sym('r%d', [1 4], 'real');
forwardVec = sym('f%d', [1 4], 'real');

indexVec = 1:8;

% ssTemp = repmat(sym('ss%d', [1 4], 'real'),4,1); 
% set which states are concentration-dependent
syms c real
sensingEdges = [2,1 ; 3, 4];

% define basic 4x4 symbolic transition matrix

RSymSimple = RFun4State(forwardVec,reverseVec);

% incorporate sensing edges      
sensingIndices = sub2ind(size(RSymSimple),sensingEdges(:,1),sensingEdges(:,2));
RSym = RSymSimple;
RSym(sensingIndices) = c*RSym(sensingIndices);

% add diagonal terms
RSym(eye(size(RSym))==1) = -sum(RSym)
%%
% %% Generate flux relations
% FluxArrayRaw = ssTemp.*RSymSimple;
% FluxArrayNet = tril(FluxArrayRaw-FluxArrayRaw');
% FluxArrayVec = FluxArrayNet(:);
% FluxEqVec = FluxArrayVec(FluxArrayVec~=0)==0;
% 
% FluxSol1 = solve(FluxEqVec
%% %%%%%%%%% derive expressions for production rate and sharpness %%%%%%%%%%

[V,D] = eig(RSym);
DLog = logical(D==0);
ssInd = find(all(DLog));

ssVecSym = V(:,ssInd) / sum(V(:,ssInd));

% full steady state vector 
ssVecFun = matlabFunction(ssVecSym,'File',[savePath 'steadyStateVecFunction4State']);

% production rate
productionRateSym = sum(ssVecSym(activeStates));
productionRateFun = matlabFunction(productionRateSym,'File',[savePath 'productionRateFunction4State'],'Optimize',false);

% take derivative wrpt c to get sharpness
sharpnessSym = diff(productionRateSym,c);
sharpnessFun = matlabFunction(sharpnessSym,'File',[savePath 'sharpnessFunction4State']);

% derive expression of variance (and, along the way, first passage times)

% generate and solve systems of equations for expected first passage time
% to each active state
stateIndex = 1:size(RSym,2);
etInfo = struct;
for CurrentState = 1:size(RSym,2)
  
    etInfo(CurrentState).CurrentState = CurrentState;
    % create symbolic vector of passage times
    etInfo(CurrentState).ETVec = sym(['ET%d' num2str(etInfo(CurrentState).CurrentState)], [1 4]);
    etInfo(CurrentState).ETVec(etInfo(CurrentState).CurrentState) = 0;
    
    % create adjusted transition matrix
    etInfo(CurrentState).RSymTo = RSym;
    etInfo(CurrentState).Rdiag = -reshape(diag(etInfo(CurrentState).RSymTo),1,[]);
    etInfo(CurrentState).RSymTo = etInfo(CurrentState).RSymTo ./ etInfo(CurrentState).Rdiag;
    etInfo(CurrentState).RSymTo(:,etInfo(CurrentState).CurrentState) = 0;
    etInfo(CurrentState).RSymTo(eye(size(etInfo(CurrentState).RSymTo))==1) = 0;
    
    % create system of equations
    etInfo(CurrentState).eqSys = etInfo(CurrentState).ETVec * etInfo(CurrentState).RSymTo + 1./etInfo(CurrentState).Rdiag;
    
    stateFilter = stateIndex~=etInfo(CurrentState).CurrentState;
    stateIndices = find(stateFilter);
    
    etInfo(CurrentState).eqSys = etInfo(CurrentState).eqSys(stateFilter);
    etInfo(CurrentState).eqSys = etInfo(CurrentState).eqSys == etInfo(CurrentState).ETVec(stateFilter);
    
    % solve
    etInfo(CurrentState).eqSol = solve(etInfo(CurrentState).eqSys,etInfo(CurrentState).ETVec(stateFilter));
    
    % get mean first passage time    
    solFields = fieldnames(etInfo(CurrentState).eqSol);
%     etInfo(CurrentState).etSolVec = sym([1 4]);
    iter = 1;
    for s = 1:length(stateIndex)
        if s == CurrentState
            etInfo(CurrentState).etSolVec(s) = ssVecSym(s)*0;
        else
            etInfo(CurrentState).etSolVec(s) = etInfo(CurrentState).eqSol.(solFields{iter});
            iter = iter + 1;
        end
    end        
    etInfo(CurrentState).ETMean = etInfo(CurrentState).etSolVec*ssVecSym;
end
 
% %%%%%%%%%%%%%%%%%%% get expressions for variance %%%%%%%%%%%%%%%%%%%%%%%


f_vec = zeros(1,size(RSym,2)); % initiation rate for each state
f_vec(activeStates) = 1; % assume all active states produce at the same rate

% construct Z matrix 
% see eqs 28 and 29 from: "Asymptotic Formulas for Markov Processes with
%                          Applications to Simulation"
ZSym = sym(size(RSym,2));
for i = 1:size(RSym,2)
    for j = 1:size(RSym,2)
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
VarFun = matlabFunction(varSum,'File',[savePath 'intrinsicVarianceFunction4State'],'Optimize',false);

% %%%%%%%%%%%%%%%%%%%%%%%%% Cycle time calculations %%%%%%%%%%%%%%%%%%%%%%%

%%% mean time to go OFF->ON %%%

% create symbolic vector of passage times
ETVecON = sym('ET%dON', [1 4]);
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
ETONMean = solVecON*ssVecSym(offStateFilter);
ETONFun = matlabFunction(ETONMean,'File',[savePath 'TauOFFFunction4State']);

%%% mean time to go ON->OFF %%%
onStateFilter = ismember(stateIndex,activeStates);

% create symbolic vector of passage times
ETVecOFF = sym('ET%dOFF', [1 4]);
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
ETOFFMean = solVecOFF*ssVecSym(onStateFilter);
ETOFFFun = matlabFunction(ETOFFMean,'File',[savePath 'TauONFunction4State']);

%%% Cycle Time %%%
TauCycleFun = matlabFunction(ETONMean+ETOFFMean,'File',[savePath 'TauCycleFunction4State']);

%%% %%%%%%%%%%%%%%%%%%%% Energy Dissipation Expressions %%%%%%%%%%%%%%%%%%%%

%%% Entropy production %%%
% See equation 5 in: Thermodynamics of Statistical Inference by Cells

entropyRateSym = 0;
ssTemp = sym('ss%d',[1,4]);
rTemp = sym('r',[4,4]);

for i = 1:size(RSym,2)
    for j = 1:size(RSym,2)
        if i ~= j && RSym(i,j)~=0
            entropyRateSym = entropyRateSym + ssVecSym(i)*RSym(j,i)*log(RSym(j,i)/RSym(i,j));
        end
    end
end

matlabFunction(entropyRateSym,'File',[savePath 'entropyRateFunction4State'],'Optimize',false);

%%% calculate flux (only well-defined for single-loop systems)
netFluxSym = ssVecSym(1)*RSym(4,1) - ssVecSym(4)*RSym(1,4);
netFluxFun = matlabFunction(netFluxSym);
matlabFunction(netFluxSym,'File',[savePath 'netFluxFunction4State']);

%%% forward and backward fluxes
[Numerator, Denominator] = numden(productionRateSym);
forwardFluxSym = RSym(1,2)*RSym(2,3)*RSym(3,4)*RSym(4,1) / Denominator;
forwardFluxFun = matlabFunction(forwardFluxSym);
matlabFunction(forwardFluxSym,'File',[savePath 'forwardFluxFunction4State']);

backwardFluxSym = RSym(2,1)*RSym(3,2)*RSym(4,3)*RSym(1,4) / Denominator;
backwardFluxFun = matlabFunction(backwardFluxSym);
matlabFunction(backwardFluxSym,'File',[savePath 'backwardFluxFunction4State']);

