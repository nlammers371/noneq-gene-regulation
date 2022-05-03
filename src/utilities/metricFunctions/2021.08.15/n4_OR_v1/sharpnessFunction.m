function sharpnessSym = sharpnessFunction(c,k12,k14,k21,k23,k32,k34,k41,k43)
%SHARPNESSFUNCTION
%    SHARPNESSSYM = SHARPNESSFUNCTION(C,K12,K14,K21,K23,K32,K34,K41,K43)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    20-Jan-2021 20:28:21

t2 = c.^2;
t3 = k12.*k14.*k23;
t4 = k14.*k21.*k23;
t5 = k14.*k21.*k32;
t6 = k12.*k14.*k43;
t7 = k12.*k23.*k41;
t8 = k14.*k21.*k43;
t9 = k12.*k34.*k41;
t10 = k14.*k32.*k43;
t11 = k12.*k41.*k43;
t12 = k23.*k34.*k41;
t13 = k32.*k34.*k41;
t14 = k32.*k41.*k43;
t17 = c.*k12.*k23.*k34;
t20 = c.*k21.*k32.*k43;
t23 = c.*k21.*k23.*k34.*2.0;
t24 = c.*k21.*k32.*k34.*2.0;
t15 = c.*t4;
t16 = c.*t5;
t18 = c.*t8;
t19 = c.*t9;
t21 = c.*t12;
t22 = c.*t13;
t25 = k21.*k23.*k34.*t2;
t26 = k21.*k32.*k34.*t2;
t27 = t3+t6+t10+t17;
t28 = t7+t11+t14+t20;
t29 = t4+t8+t12+t23;
t30 = t5+t9+t13+t24;
t31 = 1.0./t28;
t33 = t15+t18+t21+t25;
t34 = t16+t19+t22+t26;
t32 = t31.^2;
t35 = k12.*k23.*k34.*t31;
t36 = t27.*t31;
t37 = t29.*t31;
t38 = t30.*t31;
t41 = t31.*t33;
t42 = t31.*t34;
t39 = k21.*k32.*k43.*t27.*t32;
t43 = k21.*k32.*k43.*t32.*t33;
t44 = k21.*k32.*k43.*t32.*t34;
t47 = t36+t41+t42+1.0;
t40 = -t39;
t45 = -t43;
t46 = -t44;
t48 = 1.0./t47;
t49 = t48.^2;
t50 = t35+t37+t38+t40+t45+t46;
sharpnessSym = t38.*t48+t46.*t48-t49.*t50-t42.*t49.*t50;
