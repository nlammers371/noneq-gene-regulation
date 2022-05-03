function cHalf = fourStateHalfMaxGeneral(rate_vec)

for k = 1:4
    eval(['k' num2str(k) ' = rate_vec(k);'])
end
for r = 1:4
    eval(['r' num2str(r) ' = rate_vec(r+4);'])
end


cHalf = NaN;

cHalf1 = 0.5E0.*k2.^(-1).*(k3+(-0.1E1).*r1).^(-1).*r2.^(-1).*((-0.1E1).* ...
      k1.*k2.*k3+k1.*k2.*k4+(-0.1E1).*k2.*k3.*k4+k1.*k2.*r1+(-0.1E1).* ...
      k3.*r2.*r3+r1.*r2.*r3+r1.*r2.*r4+(-0.1E1).*r2.*r3.*r4+(-0.1E1).*(( ...
      -0.4E1).*k2.*(k3+(-0.1E1).*r1).*r2.*(k3.*k4.*r3+(k4+r1).*r3.*r4+ ...
      k1.*((-0.1E1).*k3.*k4+(-0.1E1).*k4.*r4+(-0.1E1).*r1.*r4))+(k2.* ...
      k3.*k4+k1.*k2.*(k3+(-0.1E1).*k4+(-0.1E1).*r1)+r2.*(k3.*r3+(-0.1E1) ...
      .*r1.*r3+(-0.1E1).*r1.*r4+r3.*r4)).^2).^(1/2));

if isreal(cHalf1)&&cHalf1>0&&~isinf(cHalf1)
    cHalf = cHalf1;
else
    cHalf2 = 0.5E0.*k2.^(-1).*(k3+(-0.1E1).*r1).^(-1).*r2.^(-1).*((-0.1E1).* ...
      k1.*k2.*k3+k1.*k2.*k4+(-0.1E1).*k2.*k3.*k4+k1.*k2.*r1+(-0.1E1).* ...
      k3.*r2.*r3+r1.*r2.*r3+r1.*r2.*r4+(-0.1E1).*r2.*r3.*r4+((-0.4E1).* ...
      k2.*(k3+(-0.1E1).*r1).*r2.*(k3.*k4.*r3+(k4+r1).*r3.*r4+k1.*(( ...
      -0.1E1).*k3.*k4+(-0.1E1).*k4.*r4+(-0.1E1).*r1.*r4))+(k2.*k3.*k4+ ...
      k1.*k2.*(k3+(-0.1E1).*k4+(-0.1E1).*r1)+r2.*(k3.*r3+(-0.1E1).*r1.* ...
      r3+(-0.1E1).*r1.*r4+r3.*r4)).^2).^(1/2));

    if isreal(cHalf2)&&cHalf2>0&&~isinf(cHalf2)
        cHalf = cHalf2;
    end
end


pd_rate = round(fourStateProductionGeneral(rate_vec,cHalf),2); 
if ~isnan(cHalf)      
    if pd_rate~=0.50
       options = optimoptions('lsqnonlin','Display','off');
       ob_fun = @(c) (fourStateProductionGeneral(rate_vec,c) - 0.5)^2;
       cHalf = lsqnonlin(ob_fun,cHalf,0,[],options);
       if round(cHalf,2) ~= 0.50
           cHalf = NaN;
       end
    end
end    



% test = ...
%   (k1.*k2.*k3+k2.*k3.*r2+r2.*r3.*(k3+r4)).*((r1.*r2+k3.*(k4+r2)).*( ...
%   k2+r3)+(r1.*r2+(k4+r1+r2).*r3).*r4+k1.*(k3.*k4+k2.*(k3+k4+r1)+(k4+ ...
%   r1).*r4)).^(-1)


% (k1.*k2.*(k4+r3)+(k2+r1).*r3.*r4).*(k2.*k3.*k4+k4.*r1.*(k3+r2)+ ...
%   r1.*r2.*r3+k1.*(k4.*(k3+r2)+r2.*r3+k2.*(k3+k4+r3))+r1.*(k3+r2).* ...
%  r4+(r1+r2).*r3.*r4+k2.*(k3+r3).*r4).^(-1);
%    
 