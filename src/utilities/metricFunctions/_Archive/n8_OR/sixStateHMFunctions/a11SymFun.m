function a11Sym = a11SymFun(b,c,cw,k12,k14,k21,k23,k32,k34,k41,k43)
%A11SYMFUN
%    A11SYM = A11SYMFUN(B,C,CW,K12,K14,K21,K23,K32,K34,K41,K43)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    21-Jan-2021 10:32:59

t2 = b.^2;
t3 = c.^2;
t4 = cw.^2;
t5 = k12.^2;
t6 = k23.^2;
t7 = k32.^2;
t8 = k43.^2;
t9 = cw.*k12.*k14.*k21.*k23.*k32;
t10 = cw.*k12.*k14.*k21.*k23.*k43;
t11 = cw.*k12.*k14.*k21.*k32.*k43;
t12 = cw.*k14.*k21.*k23.*k32.*k43;
t13 = cw.*k12.*k23.*k32.*k34.*k41;
t14 = cw.*k12.*k23.*k34.*k41.*k43;
t15 = cw.*k12.*k32.*k34.*k41.*k43;
t16 = cw.*k23.*k32.*k34.*k41.*k43;
t17 = b.*c.*k12.*k14.*k21.*k23.*k32;
t18 = b.*c.*k12.*k14.*k21.*k23.*k43;
t19 = b.*c.*k12.*k21.*k23.*k32.*k43;
t20 = b.*c.*k14.*k21.*k23.*k32.*k43;
t21 = b.*c.*k12.*k23.*k32.*k34.*k41;
t22 = b.*c.*k12.*k23.*k32.*k34.*k43;
t23 = b.*c.*k12.*k32.*k34.*k41.*k43;
t24 = b.*c.*k23.*k32.*k34.*k41.*k43;
t26 = b.*cw.*k12.*k21.*k23.*k32.*k43;
t27 = b.*cw.*k12.*k23.*k32.*k34.*k43;
t25 = b.*t10;
t28 = b.*t15;
t29 = cw.*k12.*k14.*k21.*t6;
t30 = cw.*k12.*k34.*k41.*t6;
t31 = cw.*k14.*k21.*k43.*t7;
t32 = cw.*k34.*k41.*k43.*t7;
t33 = b.*c.*k12.*k14.*k21.*t6;
t34 = b.*c.*k12.*k34.*k41.*t6;
t35 = b.*c.*k14.*k21.*k32.*t8;
t36 = b.*c.*k14.*k21.*k43.*t7;
t37 = b.*c.*k23.*k34.*k41.*t5;
t38 = b.*c.*k34.*k41.*k43.*t7;
t39 = b.*cw.*k12.*k14.*k21.*t8;
t40 = b.*cw.*k12.*k21.*k32.*t8;
t41 = b.*cw.*k14.*k21.*k32.*t8;
t42 = b.*cw.*k23.*k34.*k41.*t5;
t43 = b.*cw.*k23.*k34.*k43.*t5;
t44 = b.*cw.*k34.*k41.*k43.*t5;
t45 = b.*c.*k34.*t5.*t6;
t46 = b.*c.*k21.*t7.*t8;
t47 = b.*cw.*k34.*t5.*t6;
t48 = b.*cw.*k21.*t7.*t8;
t49 = c.*k12.*k14.*k21.*k23.*k43.*t2;
t50 = c.*k12.*k14.*k21.*k32.*k43.*t2;
t51 = c.*k12.*k23.*k34.*k41.*k43.*t2;
t52 = c.*k12.*k32.*k34.*k41.*k43.*t2;
t53 = c.*k12.*k14.*k21.*t2.*t8;
t54 = c.*k12.*k21.*k32.*t2.*t8;
t55 = c.*k23.*k34.*k43.*t2.*t5;
t56 = c.*k34.*k41.*k43.*t2.*t5;
a11Sym = (t9-t10+t11-t12+t13-t14+t15-t16+t17-t18+t19-t20+t21-t22+t23-t24-t25+t26-t27+t28-t29-t30+t31+t32-t33-t34-t35+t36+t37+t38-t39+t40-t41+t42-t43+t44-t45+t46-t47+t48-t49+t50-t51+t52-t53+t54-t55+t56+sqrt((t9-t10+t11-t12+t13-t14+t15-t16+t17-t18+t19-t20+t21-t22+t23-t24-t25+t26-t27+t28-t29-t30+t31+t32-t33-t34-t35+t36+t37+t38-t39+t40-t41+t42-t43+t44-t45+t46-t47+t48-t49+t50-t51+t52-t53+t54-t55+t56).^2-(b.*k14.*t5.*t6+b.*k14.*t7.*t8-b.*k41.*t5.*t6-b.*k41.*t7.*t8+k14.*t2.*t5.*t8-k41.*t2.*t5.*t8+b.*k12.*k14.*k32.*t8+b.*k14.*k23.*k43.*t5-b.*k12.*k32.*k41.*t8-b.*k23.*k41.*k43.*t5+k12.*k14.*k32.*t2.*t8+k14.*k23.*k43.*t2.*t5-k12.*k32.*k41.*t2.*t8-k23.*k41.*k43.*t2.*t5+b.*k12.*k14.*k23.*k32.*k43.*2.0-b.*k12.*k23.*k32.*k41.*k43.*2.0).*(k12.*k21.*k34.*t4.*t6.*4.0-k21.*k34.*k43.*t4.*t7.*4.0+c.*cw.*k12.*k21.*k34.*t6.*4.0-c.*cw.*k21.*k34.*k43.*t7.*4.0+b.*k12.*k21.*k34.*t3.*t6.*4.0-b.*k21.*k34.*k43.*t3.*t7.*4.0-k12.*k21.*k23.*k32.*k34.*t4.*4.0+k12.*k21.*k23.*k34.*k43.*t4.*4.0-k12.*k21.*k32.*k34.*k43.*t4.*4.0+k21.*k23.*k32.*k34.*k43.*t4.*4.0+b.*c.*cw.*k12.*k21.*k34.*t6.*4.0-b.*c.*cw.*k21.*k34.*k43.*t7.*4.0-c.*cw.*k12.*k21.*k23.*k32.*k34.*4.0+c.*cw.*k21.*k23.*k32.*k34.*k43.*4.0-b.*k12.*k21.*k23.*k32.*k34.*t3.*4.0+b.*k21.*k23.*k32.*k34.*k43.*t3.*4.0+k12.*k21.*k23.*k34.*k43.*t2.*t3.*4.0-k12.*k21.*k32.*k34.*k43.*t2.*t3.*4.0-b.*c.*cw.*k12.*k21.*k23.*k32.*k34.*4.0+b.*c.*cw.*k12.*k21.*k23.*k34.*k43.*8.0-b.*c.*cw.*k12.*k21.*k32.*k34.*k43.*8.0+b.*c.*cw.*k21.*k23.*k32.*k34.*k43.*4.0)))./(k12.*k21.*k34.*t4.*t6.*2.0-k21.*k34.*k43.*t4.*t7.*2.0+c.*cw.*k12.*k21.*k34.*t6.*2.0-c.*cw.*k21.*k34.*k43.*t7.*2.0+b.*k12.*k21.*k34.*t3.*t6.*2.0-b.*k21.*k34.*k43.*t3.*t7.*2.0-k12.*k21.*k23.*k32.*k34.*t4.*2.0+k12.*k21.*k23.*k34.*k43.*t4.*2.0-k12.*k21.*k32.*k34.*k43.*t4.*2.0+k21.*k23.*k32.*k34.*k43.*t4.*2.0+b.*c.*cw.*k12.*k21.*k34.*t6.*2.0-b.*c.*cw.*k21.*k34.*k43.*t7.*2.0-c.*cw.*k12.*k21.*k23.*k32.*k34.*2.0+c.*cw.*k21.*k23.*k32.*k34.*k43.*2.0-b.*k12.*k21.*k23.*k32.*k34.*t3.*2.0+b.*k21.*k23.*k32.*k34.*k43.*t3.*2.0+k12.*k21.*k23.*k34.*k43.*t2.*t3.*2.0-k12.*k21.*k32.*k34.*k43.*t2.*t3.*2.0-b.*c.*cw.*k12.*k21.*k23.*k32.*k34.*2.0+b.*c.*cw.*k12.*k21.*k23.*k34.*k43.*4.0-b.*c.*cw.*k12.*k21.*k32.*k34.*k43.*4.0+b.*c.*cw.*k21.*k23.*k32.*k34.*k43.*2.0);
