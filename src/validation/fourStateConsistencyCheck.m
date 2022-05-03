% run some basic consistency checks between 4 and 6 state functions
clear 
close all

addpath(genpath('../utilities/'))

nStates = 4;
% make sure we're pointing at the right function subfolder
rmpath(genpath('../utilities/metricFunctions/'));
addpath(['../utilities/metricFunctions/n' num2str(nStates) '/']);

c_val = 2;
rate_array_new = rand(100,8);
    
[~,newToOld] = sort([2 5 8 1 7 4 3 6]);
%           (k12,k14,k21,k23,k32,k34,k41,k43)
rate_array_old = rate_array_new;
rate_array_old = rate_array_old(:,newToOld);

% generate input cell 
valMatNew = [c_val*ones(size(rate_array_new,1),1) rate_array_new];
valCellNew = mat2cell(valMatNew,size(valMatNew,1),ones(1,size(valMatNew,2)));

%%

% calculate production rate  
tic
ProductionRateOld = fourStateProduction_v4(rate_array_old,c_val);
toc
tic
ProductionRateNew = productionRateFunction(valCellNew{:});
toc
%%
% calculate sharpness     
tic
Sharpness4 = sharpnessFunction4State(valCellOld{:});
toc
tic
Sharpness6 = sharpnessFunction(valCellNew{:});
toc
%%
% variance
tic
VarPoint4 = intrinsicVarianceFunction4State(valCellOld{:});
toc
tic
VarPoint6 = intrinsicVarianceFunction(valCellNew{:});
toc
%% 
rate_fig = figure;
scatter(ProductionRateOld,ProductionRateNew)
xlabel('production rate (4 state)')
ylabel('production rate (6 state)')
grid on

sharpnes_fig = figure;
scatter(Sharpness4,Sharpness6)
xlabel('sharpness (4 state)')
ylabel('sharpness (6 state)')
grid on

var_fig = figure;
scatter(VarPoint4,VarPoint6)
xlabel('variance (4 state)')
ylabel('variance (6 state)')
grid on
