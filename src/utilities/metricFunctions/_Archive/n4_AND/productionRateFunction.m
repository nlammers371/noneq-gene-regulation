function productionRateSym = productionRateFunction(c,k12,k14,k21,k23,k32,k34,k41,k43)
%PRODUCTIONRATEFUNCTION
%    PRODUCTIONRATESYM = PRODUCTIONRATEFUNCTION(C,K12,K14,K21,K23,K32,K34,K41,K43)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    02-Feb-2021 09:58:03

t2 = c.^2;
t3 = k12.*k23.*k41;
t4 = k12.*k41.*k43;
t5 = k32.*k41.*k43;
t6 = c.*k14.*k21.*k32;
t7 = c.*k12.*k34.*k41;
t8 = c.*k21.*k32.*k43;
t9 = c.*k32.*k34.*k41;
t10 = k21.*k32.*k34.*t2;
t11 = t3+t4+t5+t8;
t12 = 1.0./t11;
t13 = t6+t7+t9+t10;
productionRateSym = (t12.*t13)./(t12.*(k12.*k14.*k23+k12.*k14.*k43+k14.*k32.*k43+c.*k12.*k23.*k34)+t12.*t13+t12.*(c.*k14.*k21.*k23+c.*k14.*k21.*k43+c.*k23.*k34.*k41+k21.*k23.*k34.*t2)+1.0);
