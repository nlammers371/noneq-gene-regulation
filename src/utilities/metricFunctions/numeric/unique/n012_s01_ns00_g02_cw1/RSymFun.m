function RSym = RSymFun(cr,cw,a,km,kp,ka1,ki1,ka2,ki2,wa21,wi21,wa12,wi12,wpa11,wma11,wap11,wip11,wpa21,wma21,wap12,wip12)
%RSYMFUN
%    RSYM = RSYMFUN(CR,CW,A,KM,KP,KA1,KI1,KA2,KI2,WA21,WI21,WA12,WI12,WPA11,WMA11,WAP11,WIP11,WPA21,WMA21,WAP12,WIP12)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    30-Sep-2021 12:32:46

t2 = a.*km;
t3 = cr.*kp;
t4 = cw.*kp;
t5 = ka2.*wa12;
t6 = ka1.*wa21;
t7 = ka1.*wap11;
t8 = ka2.*wap12;
t9 = ki2.*wi12;
t10 = ki1.*wi21;
t11 = ki1.*wip11;
t12 = ki2.*wip12;
t13 = km.*wma11;
t14 = km.*wma21;
t15 = t2.*wma11;
t16 = t2.*wma21;
t17 = t3.*wpa11;
t18 = t3.*wpa21;
t19 = t4.*wpa11;
t20 = t4.*wpa21;
t21 = t5.*wap12;
t22 = t6.*wap11;
t23 = t9.*wip12;
t24 = t10.*wip11;
t25 = t13.*wma21;
t29 = -t7;
t30 = -t8;
t31 = -t11;
t32 = -t12;
t26 = t15.*wma21;
t27 = t17.*wpa21;
t28 = t19.*wpa21;
t33 = -t21;
t34 = -t22;
t35 = -t23;
t36 = -t24;
RSym = reshape([-ka1-ka2-t3-t4,t3,0.0,ka1,0.0,t4,ka2,0.0,0.0,0.0,0.0,0.0,km,-km+t29+t30,t7,0.0,0.0,0.0,0.0,t8,0.0,0.0,0.0,0.0,0.0,t11,-t13+t31+t33,t13,0.0,0.0,0.0,0.0,t21,0.0,0.0,0.0,ki1,0.0,t17,-ki1-t5-t17-t19,t19,0.0,0.0,0.0,0.0,t5,0.0,0.0,0.0,0.0,0.0,t15,-t15+t31+t33,t11,0.0,0.0,0.0,0.0,t21,0.0,t2,0.0,0.0,0.0,t7,-t2+t29+t30,0.0,0.0,0.0,0.0,0.0,t8,ki2,0.0,0.0,0.0,0.0,0.0,-ki2-t6-t18-t20,t18,0.0,t6,0.0,t20,0.0,t12,0.0,0.0,0.0,0.0,t14,-t14+t32+t34,t22,0.0,0.0,0.0,0.0,0.0,t23,0.0,0.0,0.0,0.0,t24,-t25+t35+t36,t25,0.0,0.0,0.0,0.0,0.0,t9,0.0,0.0,t10,0.0,t27,-t9-t10-t27-t28,t28,0.0,0.0,0.0,0.0,0.0,t23,0.0,0.0,0.0,0.0,t26,-t26+t35+t36,t24,0.0,0.0,0.0,0.0,0.0,t12,t16,0.0,0.0,0.0,t22,-t16+t32+t34],[12,12]);
