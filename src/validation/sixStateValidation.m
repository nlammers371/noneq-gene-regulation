% run some basic consistency checks between 4 and 6 state functions
clear 
close all
addpath(genpath('../utilities/'))

c_val = 1.1;
rate_array = rand(100,14);

zero_flags = [1 1 1 0 1 1 1 1 1 1 0 0 0 0 0];
%              k12,k14,k16,k21,k23,k32,k34,k41,k43,k45,k54,k56,k61,k65
inf_flags = [0 0 0 1 0 0 0 0 0 0 1 0 0 0 0];        
% [k12,k14,k21,k23,k32,k34,k41,k43]
[~,fourTosix] = sort([5 4 1 6 2 7 8 3]);

rate_array4 = rate_array(:,zero_flags(2:end)==1);
rate_array4 = rate_array4(:,fourTosix);

% generate input cell 
valMat6 = [c_val*ones(size(rate_array,1),1) rate_array];
valMat6(:,zero_flags==0) = 1e-70;
valMat6(:,inf_flags==1) = 1e70;
valCell6 = mat2cell(valMat6,size(valMat6,1),ones(1,size(valMat6,2)));

valMat4 = [c_val*ones(size(rate_array,1),1) rate_array4];
valCell4 = mat2cell(valMat4,size(valMat4,1),ones(1,size(valMat4,2)));


% calculate production rate  
tic
ProductionRate4 = productionRateFunction4State(valCell4{:});
toc
tic
ProductionRate6 = productionRateFunction(valCell6{:});
toc

% calculate sharpness     
tic
Sharpness4 = sharpnessFunction4State(valCell4{:});
toc
tic
Sharpness6 = sharpnessFunction(valCell6{:});
toc
%%
% variance
tic
VarPoint4 = intrinsicVarianceFunction4State(valCell4{:});
toc
tic
VarPoint6 = intrinsicVarianceFunction(valCell6{:});
toc
%% 
rate_fig = figure;
scatter(ProductionRate4,ProductionRate6)
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
