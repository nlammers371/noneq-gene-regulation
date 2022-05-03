function piVec = fourStateOccupancy(rate_array,c_range)

for k = 1:4
    eval(['k' num2str(k) ' = rate_array(:,k);'])
end
for r = 1:4
    eval(['r' num2str(r) ' = rate_array(:,r+4);'])
end

piVec = NaN(size(rate_array,1),4);

r2 = c_range.*r2;
k2 = c_range.*k2;

piVec(:,1) = (k1.*k2.*(k4+r1)+r1.*r2.*(k2+r3)).*((r1.*r2+k3.*(k4+r2)).*(k2+r3)+ ...
  (r1.*r2+(k4+r1+r2).*r3).*r4+k1.*(k3.*k4+k2.*(k3+k4+r1)+(k4+r1).* ...
  r4)).^(-1);

piVec(:,2) = (k1.*k3.*k4+k1.*(k4+r1).*r4+r1.*r2.*r4).*((r1.*r2+k3.*(k4+r2)).*( ...
  k2+r3)+(r1.*r2+(k4+r1+r2).*r3).*r4+k1.*(k3.*k4+k2.*(k3+k4+r1)+(k4+ ...
  r1).*r4)).^(-1);

piVec(:,3) = (k3.*k4.*(k2+r3)+(k4+r1).*r3.*r4).*((r1.*r2+k3.*(k4+r2)).*(k2+r3)+ ...
  (r1.*r2+(k4+r1+r2).*r3).*r4+k1.*(k3.*k4+k2.*(k3+k4+r1)+(k4+r1).* ...
  r4)).^(-1);

piVec(:,4) = (k1.*k2.*k3+k2.*k3.*r2+r2.*r3.*(k3+r4)).*((r1.*r2+k3.*(k4+r2)).*( ...
  k2+r3)+(r1.*r2+(k4+r1+r2).*r3).*r4+k1.*(k3.*k4+k2.*(k3+k4+r1)+(k4+ ...
  r1).*r4)).^(-1);