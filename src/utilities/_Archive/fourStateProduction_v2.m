function productionRate = fourStateProduction_v2(rate_vec,c_range,repressor_flag)

for k = 1:4
    eval(['k' num2str(k) ' = rate_vec(k);'])
end
for r = 1:4
    eval(['r' num2str(r) ' = rate_vec(r+4);'])
end
if ~repressor_flag
    k2 = k2.*c_range;
    r2 = r2.*c_range;
else
    k4 = k4.*c_range;
    r4 = r4.*c_range;
end


productionRate = ...
(k1.*k2.*k3+k2.*k3.*r2+r2.*r3.*(k3+r4)).*((r1.*r2+k3.*(k4+r2)).*( ...
  k2+r3)+(r1.*r2+(k4+r1+r2).*r3).*r4+k1.*(k3.*k4+k2.*(k3+k4+r1)+(k4+ ...
  r1).*r4)).^(-1);


% (k1.*k2.*(k4+r3)+(k2+r1).*r3.*r4).*(k2.*k3.*k4+k4.*r1.*(k3+r2)+ ...
%   r1.*r2.*r3+k1.*(k4.*(k3+r2)+r2.*r3+k2.*(k3+k4+r3))+r1.*(k3+r2).* ...
%  r4+(r1+r2).*r3.*r4+k2.*(k3+r3).*r4).^(-1);
%    
 