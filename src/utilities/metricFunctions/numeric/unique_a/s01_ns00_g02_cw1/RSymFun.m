function RSym = RSymFun(cr,cw,a,km,kp,ka,ki,wa,wi,wpa1,wma1,wap1,wip1)
%RSYMFUN
%    RSYM = RSYMFUN(CR,CW,A,KM,KP,KA,KI,WA,WI,WPA1,WMA1,WAP1,WIP1)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    30-Sep-2021 14:25:40

t2 = a.*km;
t3 = cr.*kp;
t4 = cw.*kp;
t5 = ka.*wa;
t6 = ki.*wip1;
t7 = km.*wma1;
t8 = ka.*2.0;
t9 = wma1.^2;
t10 = wpa1.^2;
t12 = ki.*wi.*2.0;
t17 = ka.*wap1.*-2.0;
t11 = t8.*wap1;
t13 = t2.*wma1;
t14 = t3.*wpa1;
t15 = t4.*wpa1;
t16 = t5.*wap1;
t18 = -t6;
t19 = t7.*wma1;
t20 = t6.*wi.*2.0;
t21 = t2.*t9;
t22 = t3.*t10;
t23 = t4.*t10;
t24 = -t16;
t25 = -t20;
RSym = reshape([-t3-t4-t8,t3,t4,t8,0.0,0.0,0.0,0.0,0.0,km,-km+t17,0.0,0.0,t11,0.0,0.0,0.0,0.0,t2,0.0,-t2+t17,0.0,0.0,t11,0.0,0.0,0.0,ki,0.0,0.0,-ki-t5-t14-t15,t14,t15,t5,0.0,0.0,0.0,t6,0.0,t7,-t7+t18+t24,0.0,0.0,t16,0.0,0.0,0.0,t6,t13,0.0,-t13+t18+t24,0.0,0.0,t16,0.0,0.0,0.0,t12,0.0,0.0,-t12-t22-t23,t22,t23,0.0,0.0,0.0,0.0,t20,0.0,t19,-t19+t25,0.0,0.0,0.0,0.0,0.0,0.0,t20,t21,0.0,-t21+t25],[9,9]);
