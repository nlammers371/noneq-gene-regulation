function RSym = RSymFun(cr,ki,ka,kp,km,wap,wip,wma,wpa,wa,wi)
%RSYMFUN
%    RSYM = RSYMFUN(CR,KI,KA,KP,KM,WAP,WIP,WMA,WPA,WA,WI)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    18-Oct-2021 13:42:15

t2 = cr.*kp;
t3 = ki.*wip;
t4 = km.*wma;
t5 = ka.*5.0;
t6 = wa.^2;
t7 = wa.^3;
t9 = wi.^2;
t10 = wi.^3;
t12 = wma.^2;
t13 = wma.^3;
t15 = wma.^5;
t16 = wpa.^2;
t17 = wpa.^3;
t19 = wpa.^5;
t20 = ka.*wa.*4.0;
t22 = ki.*wi.*2.0;
t8 = t6.^2;
t11 = t9.^2;
t14 = t12.^2;
t18 = t16.^2;
t21 = t5.*wap;
t23 = t2.*wpa;
t25 = t4.*wma;
t26 = t4.*t12;
t27 = t4.*t13;
t29 = t20.*wap;
t30 = t3.*wi.*2.0;
t31 = t2.*t16;
t32 = t2.*t17;
t34 = t2.*t19;
t36 = ka.*t6.*3.0;
t37 = ka.*t7.*2.0;
t38 = ki.*t9.*3.0;
t39 = ki.*t10.*4.0;
t43 = t3.*t9.*3.0;
t44 = t3.*t10.*4.0;
t24 = ka.*t8;
t28 = t4.*t14;
t33 = t2.*t18;
t40 = ki.*t11.*5.0;
t41 = t36.*wap;
t42 = t37.*wap;
t45 = t3.*t11.*5.0;
t35 = t24.*wap;
RSym = reshape([-t2-t5,t2,t5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,km,-km-ka.*wap.*5.0,0.0,t21,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,ki,0.0,-ki-t20-t23,t23,t20,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t3,t4,-t3-t4-ka.*wa.*wap.*4.0,0.0,t29,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t22,0.0,-t22-t31-t36,t31,t36,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t30,t25,-t25-t30-ka.*t6.*wap.*3.0,0.0,t41,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t38,0.0,-t32-t37-t38,t32,t37,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t43,t26,-t26-t43-ka.*t7.*wap.*2.0,0.0,t42,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t39,0.0,-t24-t33-t39,t33,t24,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t44,t27,-t27-t35-t44,0.0,t35,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t40,0.0,-t34-t40,t34,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t45,t28,-t28-t45],[12,12]);
