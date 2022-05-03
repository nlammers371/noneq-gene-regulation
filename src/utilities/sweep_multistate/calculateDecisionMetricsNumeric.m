function [V,T,R] = calculateDecisionMetricsNumeric(r0_vec,r1_vec,v0_vec,v1_vec,minError)


  % set threshold for decision    
  K = log((1-minError)./minError);
    
  % calculate expected drift 
  V1 = .5 .* (r0_vec - r1_vec).^2 ./ v0_vec;
  V0 = -.5 .* (r0_vec - r1_vec).^2 ./ v1_vec;
  V = .5*V1 + abs(.5*V0);
  
  % calculate diffusion
  D1 = (r0_vec-r1_vec).^2 .* v1_vec ./ (2.*v0_vec.^2);
  D0 = (r0_vec-r1_vec).^2 .* v0_vec ./ (2.*v1_vec.^2);
  
  % calculate decision times 
  B1 = V1.*K./D1;
  T1 = K./(2.*V1.*sinh(B1)) .* (exp(B1) + exp(-B1) - 2);
  
  B0 = -V0.*K./D0;
  T0 = -K./(2.*V0.*sinh(B0)) .* (exp(B0) + exp(-B0) - 2);
  
  T = .5*T1 + .5*T0;
  
  % calculate average occupancy
  R = 0.5*r1_vec + 0.5*r0_vec;  