function pdRateSymSimpSimp = pdRateSymSimpFun(k12,k14,k21,k23,k32,k34,k41,k43)
%PDRATESYMSIMPFUN
%    PDRATESYMSIMPSIMP = PDRATESYMSIMPFUN(K12,K14,K21,K23,K32,K34,K41,K43)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    19-Jan-2021 16:30:37

t2 = k14.*k21;
t3 = k12+k32;
t4 = k14+k32;
t5 = k14+k34;
t6 = k14+k41;
t7 = k23+k32;
t8 = k21+k41;
t9 = k23+k41;
t10 = k23+k43;
t11 = k21.*2.0;
t12 = k34.*2.0;
t13 = k34+t5;
t14 = k21+t8;
t15 = k32.*t6;
t16 = k34.*t9;
t17 = t2.*t10;
t18 = k34.*k41.*t3;
t19 = k34.*t14;
t20 = k21.*k32.*t13;
t21 = t2+t19;
pdRateSymSimpSimp = (k32.*(k43.*(t18+t20+k12.*k21.*t5)+k12.*k23.*t21)+k43.*(t3.*t18+k21.*k32.*(k14.*k32+k12.*t12+k32.*t12)).*1.0e+2+k12.*k23.*(t20+k34.*k41.*(k12+t3)).*1.0e+2+k12.*k43.*(t18+k21.*k32.*t5).*1.0e+4)./(k12.*(k43.*(t17+t15.*(k23+t10)+k32.*t16+k21.*k32.*t10)+k23.*t2.*(k43+t7)+k32.*k43.*(k34.*k41+k23.*(k21+k34))+k21.*k23.*k34.*t7+k23.*k34.*k41.*t7+k21.*k34.*t7.*(k43+t10)).*1.0e+2+k12.^2.*(k23.*t16+t10.*(t16+k23.*k41)+k14.*k23.*t10).*1.0e+2+t7.*(k12.*(t17+k21.*k23.*k34+k34.*t8.*t10)+k32.*k43.*t21)+k12.*k43.*(k12.*(t16+k14.*t10+k41.*t10)+k43.*(t15+k21.*t4)+t7.*(t2+k34.*t8)).*1.0e+4+k32.*k43.*(t7.*t21+k43.*(t15+t4.*t11)).*1.0e+2);
