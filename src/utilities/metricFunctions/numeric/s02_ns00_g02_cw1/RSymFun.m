function RSym = RSymFun(cr,cw,a,ki,ka,kp,km,wap,wip,wma,wpa,wmp,wa,wi)
%RSYMFUN
%    RSYM = RSYMFUN(CR,CW,A,KI,KA,KP,KM,WAP,WIP,WMA,WPA,WMP,WA,WI)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    17-Feb-2022 12:26:41

t2 = a.*km;
t3 = cr.*kp;
t4 = cw.*kp;
t5 = ka.*wa;
t6 = ki.*wip;
t7 = km.*wma;
t8 = km.*wmp;
t9 = ka.*2.0;
t10 = wap.^2;
t11 = wip.^2;
t12 = wma.^2;
t13 = wpa.^2;
t17 = ki.*wi.*2.0;
t28 = ka.*wap.*-2.0;
t14 = t3.*2.0;
t15 = t4.*2.0;
t16 = t9.*wap;
t18 = t8.*2.0;
t19 = t2.*wma;
t20 = t2.*wmp;
t21 = t3.*wpa;
t22 = t4.*wpa;
t23 = t5.*wap;
t24 = t7.*wmp;
t26 = -t3;
t27 = -t4;
t29 = -t6;
t30 = t6.*wip;
t31 = t7.*wma;
t35 = t6.*wi.*2.0;
t37 = t2.*t12;
t38 = t3.*t13;
t39 = t4.*t13;
t40 = t5.*t10;
t43 = t9.*t10;
t53 = ka.*t10.*-2.0;
t25 = t19.*wmp;
t32 = t20.*2.0;
t33 = t14.*wpa;
t34 = t15.*wpa;
t36 = t24.*2.0;
t41 = t24.*wma;
t44 = -t21;
t45 = -t22;
t46 = -t23;
t47 = -t35;
t48 = t13.*t14;
t49 = t13.*t15;
t50 = t30.*wi.*2.0;
t52 = t12.*t20;
t54 = t29.*wip;
t55 = t13.*t26;
t56 = t13.*t27;
t57 = -t40;
t42 = t25.*2.0;
t51 = t36.*wma;
t58 = -t50;
t59 = t12.*t32;
RSym = reshape([-t9-t14-t15,t15,0.0,t14,0.0,0.0,t9,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t2,-t2+t26+t27+t28,t4,0.0,t3,0.0,0.0,t16,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t32,-t32+t53,0.0,0.0,0.0,0.0,0.0,t43,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,km,0.0,0.0,-km+t26+t27+t28,t4,t3,0.0,0.0,0.0,t16,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t8,0.0,t20,-t8-t20+t53,0.0,0.0,0.0,0.0,0.0,t43,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t18,0.0,-t18+t53,0.0,0.0,0.0,0.0,0.0,t43,0.0,0.0,0.0,0.0,0.0,0.0,ki,0.0,0.0,0.0,0.0,0.0,-ki-t5-t21.*2.0-t22.*2.0,t34,0.0,t33,0.0,0.0,t5,0.0,0.0,0.0,0.0,0.0,0.0,t6,0.0,0.0,0.0,0.0,t19,-t19+t29+t44+t45+t46,t22,0.0,t21,0.0,0.0,t23,0.0,0.0,0.0,0.0,0.0,0.0,t30,0.0,0.0,0.0,0.0,t42,-t42+t54+t57,0.0,0.0,0.0,0.0,0.0,t40,0.0,0.0,0.0,0.0,0.0,0.0,t6,0.0,0.0,t7,0.0,0.0,-t7+t29+t44+t45+t46,t22,t21,0.0,0.0,0.0,t23,0.0,0.0,0.0,0.0,0.0,0.0,t30,0.0,0.0,t24,0.0,t25,-t24-t25+t54+t57,0.0,0.0,0.0,0.0,0.0,t40,0.0,0.0,0.0,0.0,0.0,0.0,t30,0.0,0.0,0.0,t36,0.0,-t36+t54+t57,0.0,0.0,0.0,0.0,0.0,t40,0.0,0.0,0.0,0.0,0.0,0.0,t17,0.0,0.0,0.0,0.0,0.0,-t17-t38.*2.0-t39.*2.0,t49,0.0,t48,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t35,0.0,0.0,0.0,0.0,t37,-t37+t47+t55+t56,t39,0.0,t38,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t50,0.0,0.0,0.0,0.0,t59,t52.*-2.0+t58,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t35,0.0,0.0,t31,0.0,0.0,-t31+t47+t55+t56,t39,t38,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t50,0.0,0.0,t41,0.0,t52,-t41-t52+t58,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t50,0.0,0.0,0.0,t51,0.0,t41.*-2.0+t58],[18,18]);
