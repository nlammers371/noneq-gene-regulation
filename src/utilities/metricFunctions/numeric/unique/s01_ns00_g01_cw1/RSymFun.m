function RSym = RSymFun(cr,cw,a,km,kp,ka1,ki1,wpa11,wma11,wap11,wip11)
%RSYMFUN
%    RSYM = RSYMFUN(CR,CW,A,KM,KP,KA1,KI1,WPA11,WMA11,WAP11,WIP11)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    30-Sep-2021 14:07:26

t2 = a.*km;
t3 = cr.*kp;
t4 = cw.*kp;
t5 = ka1.*wap11;
t6 = ki1.*wip11;
t7 = km.*wma11;
t8 = t2.*wma11;
t9 = t3.*wpa11;
t10 = t4.*wpa11;
t11 = -t5;
t12 = -t6;
RSym = reshape([-ka1-t3-t4,t3,0.0,ka1,0.0,t4,km,-km+t11,t5,0.0,0.0,0.0,0.0,t6,-t7+t12,t7,0.0,0.0,ki1,0.0,t9,-ki1-t9-t10,t10,0.0,0.0,0.0,0.0,t8,-t8+t12,t6,t2,0.0,0.0,0.0,t5,-t2+t11],[6,6]);
