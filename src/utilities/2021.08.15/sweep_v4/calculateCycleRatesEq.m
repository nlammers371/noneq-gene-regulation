function      replacementRates = calculateCycleRatesEq(rate_array, c_range, T)
    for k = 1:4
        eval(['k' num2str(k) ' = rate_array(:,k);'])
    end
    for r = 1:4
        eval(['r' num2str(r) ' = rate_array(:,r+4);'])
    end
    
    replacementRates = NaN(size(rate_array));
    
    replacementRates(:,1) = k4.^(-1).*(k4+c_range.*r2).*r3.*r4.*(c_range.*k2.*((-1)+k3.*T)+r4.*((-1)+r3.*T)).^(-1);
    replacementRates(:,5) = k2.*k3.*r2.^(-1).*(k4+c_range.*r2).*(c_range.*k2.*((-1)+k3.*T)+r4.*((-1)+r3.*T)).^(-1);
    
    replacementRates(:,2) = r1.*r4.*((-1).*c_range.*k1.*(k3+r1)+c_range.*k1.*k3.*r1.*T).^(-1).*(k1+r3+(-1).*k1.*r3.*T);
    replacementRates(:,6) = k3.*k4.*(k1+r3+(-1).*k1.*r3.*T).*((-1).*c_range.*(k3+r1).*r3+c_range.*k3.*r1.*r3.*T).^(-1);
    
    replacementRates(:,3) = k2.^(-1).*r1.*r2.*(c_range.*k2+r4).*(k4.*((-1)+k1.*T)+c_range.*r2.*((-1)+r1.*T)).^(-1);
    replacementRates(:,7) = k1.*k4.*r4.^(-1).*(c_range.*k2+r4).*(k4.*((-1)+k1.*T)+c_range.*r2.*((-1)+r1.*T)).^(-1);
    
    replacementRates(:,4) = c_range.*r2.*r3.*(k3+r1+(-1).*k3.*r1.*T).*((-1).*k3.*(k1+r3)+k1.*k3.*r3.*T).^(-1);
    replacementRates(:,8) = c_range.*k1.*k2.*(k3+r1+(-1).*k3.*r1.*T).*((-1).*r1.*(k1+r3)+k1.*r1.*r3.*T).^(-1);