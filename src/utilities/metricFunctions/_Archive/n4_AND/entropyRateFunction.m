function entropyRateSym = entropyRateFunction(c,k12,k14,k21,k23,k32,k34,k41,k43)
%ENTROPYRATEFUNCTION
%    ENTROPYRATESYM = ENTROPYRATEFUNCTION(C,K12,K14,K21,K23,K32,K34,K41,K43)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    02-Feb-2021 10:09:24

t2 = c.^2;
t3 = k12.*k14.*k23;
t4 = k12.*k14.*k43;
t5 = k12.*k23.*k41;
t6 = k14.*k32.*k43;
t7 = k12.*k41.*k43;
t8 = k32.*k41.*k43;
t9 = 1.0./c;
t10 = c.*k14.*k21.*k23;
t11 = c.*k14.*k21.*k32;
t12 = c.*k12.*k23.*k34;
t13 = c.*k14.*k21.*k43;
t14 = c.*k12.*k34.*k41;
t15 = c.*k21.*k32.*k43;
t16 = c.*k23.*k34.*k41;
t17 = c.*k32.*k34.*k41;
t18 = k21.*k23.*k34.*t2;
t19 = k21.*k32.*k34.*t2;
t20 = t3+t4+t6+t12;
t21 = t5+t7+t8+t15;
t22 = 1.0./t21;
t23 = t10+t13+t16+t18;
t24 = t11+t14+t17+t19;
t25 = t20.*t22;
t26 = t22.*t23;
t27 = t22.*t24;
t28 = t25+t26+t27+1.0;
t29 = 1.0./t28;
entropyRateSym = k14.*t29.*log(k14./k41)+c.*k34.*t29.*log((c.*k34)./k43)+k12.*t26.*t29.*log((k12.*t9)./k21)+k43.*t27.*t29.*log((k43.*t9)./k34)+k23.*t27.*t29.*log(k23./k32)+k32.*t26.*t29.*log(k32./k23)+k41.*t25.*t29.*log(k41./k14)+c.*k21.*t25.*t29.*log((c.*k21)./k12);
