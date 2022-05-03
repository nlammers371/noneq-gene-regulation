function [V,T,R,r1,r0] = calculateDecisionMetrics(param_array,simInfo)

  minError = simInfo.minError;
  
  % generate input cell arrays
  valMat1 = param_array;
  valMat1(:,simInfo.cr_index) = simInfo.cr1;
  valMat0 = param_array;
  valMat0(:,simInfo.cr_index) = simInfo.cr0;  
   
  valCellC0 = mat2cell(valMat0,size(valMat0,1),ones(1,size(valMat0,2)));
  valCellC1 = mat2cell(valMat1,size(valMat1,1),ones(1,size(valMat1,2)));
  
  % set threshold for decision    
  K = log((1-minError)./minError);

  r1 = productionRateFunction(valCellC1{:});
  r0 = productionRateFunction(valCellC0{:});
  
  v1 = intrinsicVarianceFunction(valCellC1{:});
  v0 = intrinsicVarianceFunction(valCellC0{:}); 
  
  % calculate expected drift 
  V1 = .5 .* (r0 - r1).^2 ./ v0;
  V0 = -.5 .* (r0 - r1).^2 ./ v1;
  V = .5*V1 + abs(.5*V0);
  
  % calculate diffusion
  D1 = (r0-r1).^2 .* v1 ./ (2.*v0.^2);
  D0 = (r0-r1).^2 .* v0 ./ (2.*v1.^2);
  
  % calculate decision times 
  B1 = V1.*K./D1;
  T1 = K./(2.*V1.*sinh(B1)) .* (exp(B1) + exp(-B1) - 2);
  
  B0 = -V0.*K./D0;
  T0 = -K./(2.*V0.*sinh(B0)) .* (exp(B0) + exp(-B0) - 2);
  
  T = .5*T1 + .5*T0;
  
  % calculate average occupancy
  R = 0.5*r1 + 0.5*r0;
  
%   % set threshold for decision    
%   K = log((1-minError)./minError);
% 
%   r1 = fourStateProduction_v4(rate_array,c1);
%   r0 = fourStateProduction_v4(rate_array,c0);
%   
%   v1 = fourStateVariance_v4(rate_array,c1);
%   v0 = fourStateVariance_v4(rate_array,c0);