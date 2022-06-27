function out1 = TauCycleFunction(cr,kip,kap,kim,kam,kpi,kmi,kpa,kma)
%TAUCYCLEFUNCTION
%    OUT1 = TAUCYCLEFUNCTION(CR,KIP,KAP,KIM,KAM,KPI,KMI,KPA,KMA)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    25-Apr-2021 12:57:42

t2 = cr.*kpa;
t3 = cr.*kpi;
t4 = kam.*kap;
t5 = kam.*kmi;
t6 = kim.*kip;
t7 = kim.*kma;
t8 = cr.^2;
t9 = kap.*t3;
t10 = kip.*t2;
t11 = kma.*t4;
t12 = kap.*t7;
t13 = kip.*t5;
t14 = kma.*t5;
t15 = kmi.*t6;
t16 = kmi.*t7;
t17 = t2.*t4;
t20 = t2.*t5;
t22 = t3.*t6;
t23 = t3.*t7;
t18 = kam.*t10;
t19 = kim.*t9;
t21 = kma.*t9;
t24 = kmi.*t10;
t25 = t2.*t9;
t26 = t3.*t10;
t27 = t4+t5+t9;
t28 = t6+t7+t10;
t29 = 1.0./t27;
t30 = 1.0./t28;
t31 = t11+t13+t14+t21;
t32 = t12+t15+t16+t24;
t34 = t17+t19+t20+t25;
t35 = t18+t22+t23+t26;
t33 = 1.0./t31;
t36 = t32.*t33;
t37 = t33.*t34;
t38 = t33.*t35;
t39 = t36+t37+t38+1.0;
t40 = 1.0./t39;
out1 = (kam.*t30.*t36.*t40.*(kip+kma+t2)+kap.*t30.*t38.*t40.*(kim+kma+t2))./(kam.*t36.*t40+kap.*t38.*t40)+(kim.*t29.*t40.*(kap+kmi+t3)+kip.*t29.*t37.*t40.*(kam+kmi+t3))./(kim.*t40+kip.*t37.*t40);