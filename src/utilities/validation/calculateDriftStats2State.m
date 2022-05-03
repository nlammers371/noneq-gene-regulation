function [r0, r1, rt, v0, v1, vt, V_mRNA, D_mRNA, V_micro, D_micro, tau_micro, tau_mRNA]...
                          = calculateDriftStats2State(kon,koff,c0,c1,c_true,K)

  r0 = c0*kon / (c0*kon + koff);
  r1 = c1*kon / (c1*kon + koff);
  rt = c_true*kon / (c_true*kon + koff);
  
  v0 = 2*c0*kon*koff / (c0*kon+koff)^3;
  v1 = 2*c1*kon*koff / (c1*kon+koff)^3;
  vt = 2*c_true*kon*koff / (c_true*kon+koff)^3;

  % calculate expected drift and diffusion rates for mRNA-driven decision
  V_mRNA = .5*( (r0 - rt).^2 ./ v0 - (r1 - rt).^2 ./ v1); 
  D_mRNA = (((r1 - rt).*v0 + (rt-r0).*v1).^2.*vt) /(2*v1.^2 .* v0.^2);

  % calculate expected drift and diffusion rates for transition-driven decision
  V_micro = koff.*kon ./ (koff + c_true*kon).*(c0-c1+c_true.*log(c1/c0));
  
  D_micro = koff.*c_true.*kon ./ (koff + c_true*kon).^3 * ( (kon.*(c1-c0)).^2+.5*log(c1./c0).^2.*...
          (koff.^2+kon.^2.*c_true.^2)+kon.*(c1-c0).*log(c1./c0).*(koff-c_true.*kon));
        
        
  % calculate expected decision times
  B_micro = V_micro*K/D_micro;
  tau_micro = K/(2*V_micro*sinh(B_micro)) * (exp(B_micro) + exp(-B_micro) - 2);
  
  B_mRNA = V_mRNA*K/D_mRNA;
  tau_mRNA = K/(2*V_mRNA*sinh(B_mRNA)) * (exp(B_mRNA) + exp(-B_mRNA) - 2);
          
  