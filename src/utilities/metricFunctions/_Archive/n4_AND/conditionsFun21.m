function out1 = conditionsFun21(c,k12,k14,k21,k23,k32,k34,k41,k43)
%CONDITIONSFUN21
%    OUT1 = CONDITIONSFUN21(C,K12,K14,K21,K23,K32,K34,K41,K43)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    02-Feb-2021 10:03:20

t2 = c.^2;
t3 = k12.^2;
t4 = k14.^2;
t5 = k21.^2;
t6 = k23.^2;
t7 = k32.^2;
t8 = k34.^2;
t9 = k41.^2;
t10 = k43.^2;
t11 = k12.*k14.*k23;
t12 = k12.*k14.*k43;
t13 = k12.*k23.*k41;
t14 = k14.*k32.*k43;
t15 = k12.*k41.*k43;
t16 = k32.*k41.*k43;
t17 = c.*k12.*k23.*k34;
t18 = c.*k14.*k21.*k43;
t19 = c.*k12.*k34.*k41;
t20 = c.*k21.*k32.*k43;
t21 = -t19;
t22 = t11.*t16.*4.0;
t23 = t11.*t20.*2.0;
t24 = c.*k21.*k41.*k43.*t11.*2.0;
t25 = c.*k32.*k34.*k43.*t11.*2.0;
t26 = c.*k21.*k32.*k41.*t12.*4.0;
t27 = c.*k34.*k41.*k43.*t11.*4.0;
t28 = t13.*t20.*2.0;
t29 = c.*k32.*k34.*k41.*t12.*2.0;
t30 = c.*k32.*k34.*k43.*t13.*2.0;
t31 = t12+t15;
t32 = k12.*k23.*k32.*k43.*t4.*2.0;
t33 = k12.*k23.*k32.*k43.*t9.*2.0;
t34 = c.*k12.*k21.*k23.*k43.*t4.*2.0;
t35 = k12.*t4.*t20.*4.0;
t36 = c.*k14.*k23.*k34.*k41.*t3.*2.0;
t37 = c.*k14.*k21.*k32.*k41.*t10.*2.0;
t38 = k43.*t9.*t17.*4.0;
t39 = c.*k12.*k32.*k34.*k43.*t9.*2.0;
t42 = t3.*t4.*t6;
t43 = t3.*t6.*t9;
t44 = t4.*t7.*t10;
t45 = t7.*t9.*t10;
t46 = k14.*k41.*t3.*t6.*2.0;
t47 = k14.*k41.*t7.*t10.*2.0;
t48 = c.*k14.*k34.*t3.*t6.*2.0;
t49 = c.*k14.*k21.*t7.*t10.*2.0;
t50 = c.*k23.*k34.*t3.*t9.*2.0;
t51 = c.*k34.*k41.*t3.*t6.*2.0;
t52 = c.*k21.*k32.*t4.*t10.*2.0;
t53 = c.*k21.*k41.*t7.*t10.*2.0;
t57 = k21.*k34.*k43.*t2.*t11.*2.0;
t58 = k21.*k32.*k34.*t2.*t12.*4.0;
t59 = k21.*k34.*k41.*t2.*t12.*2.0;
t60 = k12.*k21.*k23.*k32.*k34.*k43.*t2.*2.0;
t61 = k21.*k34.*k43.*t2.*t13.*4.0;
t62 = k21.*k32.*k34.*t2.*t15.*2.0;
t68 = t2.*t3.*t6.*t8;
t69 = t2.*t4.*t5.*t10;
t70 = t2.*t3.*t8.*t9;
t71 = t2.*t5.*t7.*t10;
t72 = k14.*k32.*t2.*t5.*t10.*2.0;
t73 = k23.*k41.*t2.*t3.*t8.*2.0;
t75 = t13+t16+t20;
t40 = -t24;
t41 = -t27;
t54 = -t34;
t55 = -t36;
t56 = -t38;
t63 = t31.^2;
t64 = -t50;
t65 = -t57;
t66 = -t59;
t67 = -t61;
t74 = -t73;
t76 = t31.*t75.*2.0;
t77 = t22+t23+t25+t26+t28+t29+t30+t32+t33+t35+t37+t39+t40+t41+t42+t43+t44+t45+t46+t47+t48+t49+t51+t52+t53+t54+t55+t56+t58+t60+t62+t64+t65+t66+t67+t68+t69+t70+t71+t72+t74;
t78 = sqrt(t77);
t79 = -t78;
t80 = t11+t14+t17+t18+t21+t75+t79;
t81 = t80.^2;
t82 = t15.*t80;
out1 = ((((((t13.*t31.*2.0+t16.*t31.*2.0+t20.*t31.*2.0 ~= t82) & (t76 ~= t15.*(t11+t14+t17+t18+t21+t75+t78))) & (t76 ~= t82)) & (t19+t78 ~= t11+t14+t17+t18+t75)) & (t12.*t81+t15.*t81+c.*k14.*k21.*k23.*t63.*4.0+c.*k14.*k21.*k32.*t63.*4.0+c.*k23.*k34.*k41.*t63.*4.0+c.*k32.*k34.*k41.*t63.*4.0+k21.*k23.*k34.*t2.*t63.*4.0+k21.*k32.*k34.*t2.*t63.*4.0 ~= t11.*t31.*t80.*2.0+t13.*t31.*t80.*2.0+t14.*t31.*t80.*2.0+t16.*t31.*t80.*2.0+t17.*t31.*t80.*2.0+t18.*t31.*t80.*2.0+t19.*t31.*t80.*2.0+t20.*t31.*t80.*2.0)) & (t11+t14+t17+t18+t75+t78 ~= t19));
