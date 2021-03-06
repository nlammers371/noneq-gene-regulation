function ETOFFMean = TauONFunction(cr,kip,kap,kim,kam,kpi,kmi,kpa,kma)
%TAUONFUNCTION
%    ETOFFMEAN = TAUONFUNCTION(CR,KIP,KAP,KIM,KAM,KPI,KMI,KPA,KMA)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    25-Apr-2021 12:57:35

t2 = cr.*kpa;
t3 = kim.*kip;
t4 = kim.*kma;
t5 = cr.^2;
t7 = kam.*kap.*kma;
t9 = kam.*kip.*kmi;
t10 = kam.*kma.*kmi;
t15 = cr.*kap.*kim.*kpi;
t17 = cr.*kap.*kma.*kpi;
t6 = kip.*t2;
t8 = kap.*t4;
t11 = kmi.*t3;
t12 = kmi.*t4;
t13 = kam.*kap.*t2;
t16 = kam.*kmi.*t2;
t18 = cr.*kpi.*t3;
t19 = cr.*kpi.*t4;
t21 = cr.*kap.*kpi.*t2;
t25 = t7+t9+t10+t17;
t14 = kam.*t6;
t20 = kmi.*t6;
t22 = cr.*kpi.*t6;
t23 = t3+t4+t6;
t27 = 1.0./t25;
t28 = t13+t15+t16+t21;
t24 = 1.0./t23;
t26 = t8+t11+t12+t20;
t29 = t14+t18+t19+t22;
t31 = t27.*t28;
t30 = t26.*t27;
t32 = t27.*t29;
t33 = t30+t31+t32+1.0;
t34 = 1.0./t33;
ETOFFMean = (kam.*t24.*t30.*t34.*(kip+kma+t2)+kap.*t24.*t32.*t34.*(kim+kma+t2))./(kam.*t30.*t34+kap.*t32.*t34);
