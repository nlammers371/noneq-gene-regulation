function RSym = RSymFun(cr,ki,ka,kp,km,wap,wip,wma,wpa,wmp,wa,wi)
%RSYMFUN
%    RSYM = RSYMFUN(CR,KI,KA,KP,KM,WAP,WIP,WMA,WPA,WMP,WA,WI)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    18-Oct-2021 13:40:37

t2 = cr.*kp;
t3 = ki.*wip;
t4 = km.*wma;
t5 = ka.*3.0;
t6 = wa.^2;
t7 = wap.^2;
t8 = wap.^3;
t9 = wi.^2;
t10 = wip.^2;
t11 = wip.^3;
t12 = wma.^2;
t13 = wma.^3;
t14 = wmp.^2;
t15 = wpa.^2;
t16 = wpa.^3;
t19 = ka.*wa.*2.0;
t21 = ki.*wi.*2.0;
t22 = km.*wmp.*2.0;
t17 = t2.*2.0;
t18 = t2.*3.0;
t20 = t5.*wap;
t23 = t2.*wpa;
t24 = ka.*t6;
t25 = t3.*wip;
t26 = t3.*t10;
t27 = t4.*wma;
t28 = t4.*t12;
t31 = t19.*wap;
t32 = t3.*wi.*2.0;
t33 = t4.*wmp.*2.0;
t34 = t2.*t15;
t35 = t2.*t16;
t37 = t5.*t7;
t38 = t5.*t8;
t39 = ki.*t9.*3.0;
t40 = km.*t14.*3.0;
t45 = t7.*t19;
t46 = t8.*t19;
t49 = t3.*t9.*3.0;
t51 = t4.*t14.*3.0;
t29 = t17.*wpa;
t30 = t18.*wpa;
t36 = t24.*wap;
t41 = t15.*t17;
t42 = t15.*t18;
t43 = t16.*t17;
t44 = t16.*t18;
t47 = t25.*wi.*2.0;
t48 = t26.*wi.*2.0;
t50 = t27.*wmp.*2.0;
t52 = t28.*wmp.*2.0;
t53 = t7.*t24;
t54 = t8.*t24;
t55 = t9.*t25.*3.0;
t56 = t9.*t26.*3.0;
t57 = t14.*t27.*3.0;
t58 = t14.*t28.*3.0;
RSym = reshape([-t5-t18,t18,0.0,0.0,t5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,km,-km-t17-ka.*wap.*3.0,t17,0.0,0.0,t20,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t22,-t2-t22-ka.*t7.*3.0,t2,0.0,0.0,t37,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t40,-t40-ka.*t8.*3.0,0.0,0.0,0.0,t38,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,ki,0.0,0.0,0.0,-ki-t19-t23.*3.0,t30,0.0,0.0,t19,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t3,0.0,0.0,t4,-t3-t4-t23.*2.0-ka.*wa.*wap.*2.0,t29,0.0,0.0,t31,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t25,0.0,0.0,t33,-t23-t25-t33-ka.*t7.*wa.*2.0,t23,0.0,0.0,t45,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t26,0.0,0.0,t51,-t26-t51-ka.*t8.*wa.*2.0,0.0,0.0,0.0,t46,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t21,0.0,0.0,0.0,-t21-t24-t34.*3.0,t42,0.0,0.0,t24,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t32,0.0,0.0,t27,-t27-t32-t34.*2.0-t36,t41,0.0,0.0,t36,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t47,0.0,0.0,t50,-t34-t47-t50-t53,t34,0.0,0.0,t53,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t48,0.0,0.0,t57,-t48-t54-t57,0.0,0.0,0.0,t54,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t39,0.0,0.0,0.0,t35.*-3.0-t39,t44,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t49,0.0,0.0,t28,-t28-t35.*2.0-t49,t43,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t55,0.0,0.0,t52,-t35-t52-t55,t35,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t56,0.0,0.0,t58,-t56-t58],[16,16]);
