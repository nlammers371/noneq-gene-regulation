function productionRate = fourStateProduction(rate_vec,c_range,frequency_flag)

for k = 1:4
    eval(['k' num2str(k) ' = rate_vec(k);'])
end
for r = 1:4
    eval(['r' num2str(r) ' = rate_vec(r+4);'])
end
if frequency_flag
    k1 = k1.*c_range;
    r3 = r3.*c_range;
else
    k3 = k3./c_range;
    r1 = r1./c_range;
end

% % generate matrix
% R = [-k1-r4       r1             0          k4; 
%        k1       -r1-k2          r2          0
%         0         k2          -r2-k3        r3
%        r4          0            k3        -r3-k4 ];
%    
% % calculate eigenvalues and vectors
% [V,D] = eig(R);
% % find minimum eigenvale
% [~,mi] = max(real(diag(D)));
% ss = V(:,mi)/sum(V(:,mi));
% productionRate = ss(3)

productionRate = (k1.*k2.*(k4+r3)+(k2+r1).*r3.*r4).*(k2.*k3.*k4+k4.*r1.*(k3+r2)+ ...
  r1.*r2.*r3+k1.*(k4.*(k3+r2)+r2.*r3+k2.*(k3+k4+r3))+r1.*(k3+r2).* ...
 r4+(r1+r2).*r3.*r4+k2.*(k3+r3).*r4).^(-1);
   
 