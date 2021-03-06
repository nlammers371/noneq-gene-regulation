function [TauOn,TauOff,KonEff,KoffEff,dKondC,dKoffdC] = fourStateCycleTime_v4(rate_array,c_range)

for k = 1:4
        eval(['k' num2str(k) ' = rate_array(:,k);'])
end
for r = 1:4
    eval(['r' num2str(r) ' = rate_array(:,r+4);'])
end


KonEff = (c_range*k1.*k2.*k3.*(k4+r1)+c_range^2.*k2.*k3.*r1.*r2+k1.*k3.*k4.*r3+k1.*( ...
  k4+r1).*r3.*r4+c_range*r1.*r2.*r3.*(k3+r4)).*(k1.*k3.*k4+c_range*k1.*k2.*( ...
  k4+r1)+c_range^2.*k2.*r1.*r2+k1.*(k4+r1).*r4+c_range*r1.*r2.*(r3+r4)).^(-1); 

KoffEff = (c_range*k2.*k3.*(k1+k4+c_range*r2)+k3.*(k4+c_range*r2).*r3+(k4+r1+c_range*r2).*r3.* ...
  r4).^(-1).*(c_range*k1.*k2.*k3.*(k4+r1)+c_range^2.*k2.*k3.*r1.*r2+k1.*k3.* ...
  k4.*r3+k1.*(k4+r1).*r3.*r4+c_range*r1.*r2.*r3.*(k3+r4));

dKondC = (k3+(-1).*r3).*(k1.*k3.*k4.*(k1.*k2.*(k4+r1)+r1.*r2.*(2.*c_range*k2+r3) ...
  )+(k2.*(k1.*(k4+r1)+c_range*r1.*r2).^2+k1.*r1.*(k4+r1).*r2.*r3).*r4).*( ...
  k1.*k3.*k4+c_range*k1.*k2.*(k4+r1)+c_range^2.*k2.*r1.*r2+k1.*(k4+r1).*r4+c_range* ...
  r1.*r2.*(r3+r4)).^(-2);

dKoffdC = (c_range*k2.*k3.*(k1+k4+c_range*r2)+k3.*(k4+c_range*r2).*r3+(k4+r1+c_range*r2).*r3.* ...
  r4).^(-2).*(c_range^2.*k2.^2.*k3.^2.*k4.*((-1).*k1+r1).*r2+(-2).*c_range* ...
  k2.*k3.*(k1+(-1).*r1).*r2.*r3.*(k3.*k4+(k4+r1).*r4)+(-1).*(k1+(-1) ...
  .*r1).*r3.*(k3.*k4+(k4+r1).*r4).*(k1.*k2.*k3+r2.*r3.*(k3+r4)));

TauOff = (k2.*k3+r3.*(k3+r4)).^(-1).*((r1.*r2+k3.*(k4+r2)).*(k2+r3)+(r1.* ...
  r2+(k4+r1+r2).*r3).*r4+k1.*(k3.*k4+k2.*(k3+k4+r1)+(k4+r1).*r4)).^( ...
  -1).*(r1.*r2.*((k2+r3).^2+(2.*k2+k3+r3).*r4+r4.^2)+k1.*(k2.^2.*( ...
  k4+r1)+(k3+r4).*(k3.*k4+(k4+r1).*r4)+k2.*(k3.*k4+(k4+r1).*(r3+2.* ...
  r4))));

TauOn = (k1.*(k4+r1)+r1.*r2).^(-1).*(k3.*(k1.^2.*k2+(k4.*(k4+r1)+2.*k4.* ...
  r2+r2.^2).*(k2+r3)+k1.*(k2.*(k4+2.*r2)+r2.*r3))+((k4+r1).^2+(k1+ ...
  2.*k4+r1).*r2+r2.^2).*r3.*r4).*((r1.*r2+k3.*(k4+r2)).*(k2+r3)+( ...
  r1.*r2+(k4+r1+r2).*r3).*r4+k1.*(k3.*k4+k2.*(k3+k4+r1)+(k4+r1).*r4) ...
  ).^(-1);

% TauOff = ...
%   (c_range*k1.*k2.*k3.*(k4+r1)+c_range^2.*k2.*k3.*r1.*r2+k1.*k3.*k4.*r3+k1.*( ...
%   k4+r1).*r3.*r4+c_range*r1.*r2.*r3.*(k3+r4)).^(-1).*(k1.*k3.*k4+c_range*k1.* ...
%   k2.*(k4+r1)+c_range^2.*k2.*r1.*r2+k1.*(k4+r1).*r4+c_range*r1.*r2.*(r3+r4));
% 
% TauOn = (c_range*k2.*k3.*(k1+k4+c_range*r2)+k3.*(k4+c_range*r2).*r3+(k4+r1+c_range*r2).*r3.* ...
%   r4).*(c_range*k1.*k2.*k3.*(k4+r1)+c_range^2.*k2.*k3.*r1.*r2+k1.*k3.*k4.*r3+ ...
%   k1.*(k4+r1).*r3.*r4+c_range*r1.*r2.*r3.*(k3+r4)).^(-1);


