function RSym = RSymFun(cr,ki,ka,kp,km,wap,wip,wma,wpa,wmp)
%RSYMFUN
%    RSYM = RSYMFUN(CR,KI,KA,KP,KM,WAP,WIP,WMA,WPA,WMP)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    18-Oct-2021 13:40:05

t2 = cr.*kp;
t3 = ka.*wap;
t4 = ki.*wip;
t5 = km.*wma;
t6 = wap.^2;
t7 = wap.^3;
t8 = wip.^2;
t9 = wip.^3;
t10 = wmp.^2;
t13 = km.*wmp.*2.0;
t11 = t2.*2.0;
t12 = t2.*3.0;
t14 = t2.*wpa;
t15 = t3.*wap;
t16 = t3.*t6;
t17 = t4.*wip;
t18 = t4.*t8;
t21 = t5.*wmp.*2.0;
t22 = km.*t10.*3.0;
t23 = t5.*t10.*3.0;
t19 = t11.*wpa;
t20 = t12.*wpa;
RSym = reshape([-ka-t12,t12,0.0,0.0,ka,0.0,0.0,0.0,km,-km-t3-t11,t11,0.0,0.0,t3,0.0,0.0,0.0,t13,-t2-t13-t15,t2,0.0,0.0,t15,0.0,0.0,0.0,t22,-t16-t22,0.0,0.0,0.0,t16,ki,0.0,0.0,0.0,-ki-t14.*3.0,t20,0.0,0.0,0.0,t4,0.0,0.0,t5,-t4-t5-t14.*2.0,t19,0.0,0.0,0.0,t17,0.0,0.0,t21,-t14-t17-t21,t14,0.0,0.0,0.0,t18,0.0,0.0,t23,-t18-t23],[8,8]);
