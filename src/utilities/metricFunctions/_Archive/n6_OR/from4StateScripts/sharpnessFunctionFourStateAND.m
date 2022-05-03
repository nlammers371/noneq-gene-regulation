function sharpnessSym = sharpnessFunctionFourStateAND(c,k12,k14,k21,k23,k32,k34,k41,k43)
%SHARPNESSFUNCTION
%    SHARPNESSSYM = SHARPNESSFUNCTION(C,K12,K14,K21,K23,K32,K34,K41,K43)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    02-Feb-2021 09:58:04

t2 = c.^2;
t3 = k12.*k14.*k23;
t4 = k14.*k21.*k32;
t5 = k12.*k14.*k43;
t6 = k12.*k23.*k41;
t7 = k12.*k34.*k41;
t8 = k14.*k32.*k43;
t9 = k12.*k41.*k43;
t10 = k32.*k34.*k41;
t11 = k32.*k41.*k43;
t12 = c.*k14.*k21.*k23;
t14 = c.*k12.*k23.*k34;
t15 = c.*k14.*k21.*k43;
t17 = c.*k21.*k32.*k43;
t18 = c.*k23.*k34.*k41;
t20 = c.*k21.*k32.*k34.*2.0;
t13 = c.*t4;
t16 = c.*t7;
t19 = c.*t10;
t21 = k21.*k23.*k34.*t2;
t22 = k21.*k32.*k34.*t2;
t23 = t3+t5+t8+t14;
t24 = t6+t9+t11+t17;
t25 = t4+t7+t10+t20;
t26 = 1.0./t24;
t28 = t12+t15+t18+t21;
t29 = t13+t16+t19+t22;
t27 = t26.^2;
t30 = t23.*t26;
t31 = t26.*t28;
t32 = t26.*t29;
t33 = t30+t31+t32+1.0;
t34 = 1.0./t33;
sharpnessSym = -t32.*t34.^2.*(t26.*(k14.*k21.*k23+k14.*k21.*k43+k23.*k34.*k41+c.*k21.*k23.*k34.*2.0)+t25.*t26+k12.*k23.*k34.*t26-k21.*k32.*k43.*t23.*t27-k21.*k32.*k43.*t27.*t28-k21.*k32.*k43.*t27.*t29)+t25.*t26.*t34-k21.*k32.*k43.*t27.*t29.*t34;
