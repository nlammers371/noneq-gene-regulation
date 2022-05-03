function k3Half = k3HM(rate_array,c_range)

for k = 1:4
    eval(['k' num2str(k) ' = rate_array(:,k);'])
end
for r = 1:4
    eval(['r' num2str(r) ' = rate_array(:,r+4);'])
end

% k1Min = (((k4+r1+c_range.*r2).*r3.*r4) - c_range.*r1.*r2.*(c_range.*k2+r3+r4))./((k4+r1).*(c_range.*k2+r4));
% 
k3Half = (c_range.*k1.*k2+(-1).*k1.*k4+c_range.*k2.*k4+c_range.^2.*k2.*r2+k4.*r3+c_range.*r2.*r3) ...
  .^(-1).*(c_range.*k1.*k2.*k4+c_range.*k1.*k2.*r1+c_range.^2.*k2.*r1.*r2+c_range.*r1.*r2.* ...
  r3+k1.*k4.*r4+k1.*r1.*r4+c_range.*r1.*r2.*r4+(-1).*k4.*r3.*r4+(-1).*r1.* ...
  r3.*r4+(-1).*c_range.*r2.*r3.*r4);






% (k1.*k2.*(k4+r3)+(k2+r1).*r3.*r4).*(k2.*k3.*k4+k4.*r1.*(k3+r2)+ ...
%   r1.*r2.*r3+k1.*(k4.*(k3+r2)+r2.*r3+k2.*(k3+k4+r3))+r1.*(k3+r2).* ...
%  r4+(r1+r2).*r3.*r4+k2.*(k3+r3).*r4).^(-1);
%    
 