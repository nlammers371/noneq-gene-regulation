function RSym = RSymFun(cr,ki,ka,kp,km,wap,wip,wma,wpa,wmp)
%RSYMFUN
%    RSYM = RSYMFUN(CR,KI,KA,KP,KM,WAP,WIP,WMA,WPA,WMP)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    18-Oct-2021 13:40:03

t2 = cr.*kp;
t3 = ka.*wap;
t4 = ki.*wip;
t5 = km.*wma;
t6 = wap.^2;
t7 = wip.^2;
t9 = km.*wmp.*2.0;
t8 = t2.*2.0;
t10 = t2.*wpa;
t11 = t3.*wap;
t12 = t4.*wip;
t14 = t5.*wmp.*2.0;
t13 = t8.*wpa;
RSym = reshape([-ka-t8,t8,0.0,ka,0.0,0.0,km,-km-t2-t3,t2,0.0,t3,0.0,0.0,t9,-t9-t11,0.0,0.0,t11,ki,0.0,0.0,-ki-t10.*2.0,t13,0.0,0.0,t4,0.0,t5,-t4-t5-t10,t10,0.0,0.0,t12,0.0,t14,-t12-t14],[6,6]);
