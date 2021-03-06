function RSym = RSymFun(cr,cw,b,kap1,kam1,kip1,kim1,kap2,kam2,kip2,kim2,kpi,kpa,kmi,kma)
%RSYMFUN
%    RSYM = RSYMFUN(CR,CW,B,KAP1,KAM1,KIP1,KIM1,KAP2,KAM2,KIP2,KIM2,KPI,KPA,KMI,KMA)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    03-Aug-2021 16:20:28

t2 = b.*kma;
t3 = b.*kmi;
t4 = cr.*kpa;
t5 = cr.*kpi;
t6 = cw.*kpa;
t7 = cw.*kpi;
t8 = -kam1;
t9 = -kam2;
t10 = -kap1;
t11 = -kap2;
t12 = -kim1;
t13 = -kim2;
t14 = -kip1;
t15 = -kip2;
t16 = -kma;
t17 = -kmi;
t18 = -t2;
t19 = -t3;
t20 = -t4;
t21 = -t5;
t22 = -t6;
t23 = -t7;
RSym = reshape([t8+t9+t20+t22,t4,0.0,kam1,0.0,t6,kam2,0.0,0.0,0.0,0.0,0.0,kmi,t10+t11+t17,kap1,0.0,0.0,0.0,0.0,kap2,0.0,0.0,0.0,0.0,0.0,kip1,t11+t14+t17,kmi,0.0,0.0,0.0,0.0,kap2,0.0,0.0,0.0,kim1,0.0,t5,t9+t12+t21+t23,t7,0.0,0.0,0.0,0.0,kam2,0.0,0.0,0.0,0.0,0.0,t3,t11+t14+t19,kip1,0.0,0.0,0.0,0.0,kap2,0.0,t3,0.0,0.0,0.0,kap1,t10+t11+t19,0.0,0.0,0.0,0.0,0.0,kap2,kim2,0.0,0.0,0.0,0.0,0.0,t8+t13+t20+t22,t4,0.0,kam1,0.0,t6,0.0,kip2,0.0,0.0,0.0,0.0,kma,t10+t15+t16,kap1,0.0,0.0,0.0,0.0,0.0,kip2,0.0,0.0,0.0,0.0,kip1,t14+t15+t16,kma,0.0,0.0,0.0,0.0,0.0,kim2,0.0,0.0,kim1,0.0,t5,t12+t13+t21+t23,t7,0.0,0.0,0.0,0.0,0.0,kip2,0.0,0.0,0.0,0.0,t2,t14+t15+t18,kip1,0.0,0.0,0.0,0.0,0.0,kip2,t2,0.0,0.0,0.0,kap1,t10+t15+t18],[12,12]);
