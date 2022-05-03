function out1 = hmFun31(cr,kip,kap,kim,kam,kpi,kmi,kpa,kma)
%HMFUN31
%    OUT1 = HMFUN31(CR,KIP,KAP,KIM,KAM,KPI,KMI,KPA,KMA)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    26-Apr-2021 13:11:32

t2 = cr.^2;
t3 = cr.^3;
t4 = kam.^2;
t5 = kap.^2;
t6 = kim.^2;
t7 = kip.^2;
t8 = kma.^2;
t9 = kmi.^2;
t10 = kpa.^2;
t11 = kpi.^2;
out1 = ((sqrt(t5.*t6.*t8+t4.*t7.*t9+t6.*t8.*t9+kap.*kmi.*t6.*t8.*2.0+t2.*t4.*t7.*t10+t2.*t5.*t6.*t11+t2.*t6.*t8.*t11+t2.*t7.*t9.*t10+t2.^2.*t7.*t10.*t11+kam.*kim.*kip.*kma.*t9.*2.0-cr.*kam.*kpa.*t7.*t9.*2.0+cr.*kap.*kpi.*t6.*t8.*2.0-cr.*kma.*kpi.*t5.*t6.*2.0-cr.*kmi.*kpa.*t4.*t7.*2.0+cr.*kmi.*kpi.*t6.*t8.*2.0-kap.*kma.*t2.*t6.*t11.*2.0+kam.*kmi.*t2.*t7.*t10.*2.0+kam.*kpi.*t3.*t7.*t10.*2.0+kmi.*kpi.*t3.*t7.*t10.*2.0+kam.*kap.*kim.*kip.*kma.*kmi.*2.0+cr.*kam.*kim.*kip.*kpa.*t9.*4.0-cr.*kap.*kma.*kmi.*kpi.*t6.*2.0+cr.*kim.*kip.*kma.*kpa.*t9.*2.0+kap.*kim.*kip.*kma.*t2.*t11.*4.0+kap.*kim.*kip.*kpa.*t3.*t11.*2.0-kam.*kmi.*kpa.*kpi.*t2.*t7.*2.0+kim.*kip.*kma.*kpa.*t3.*t11.*2.0+cr.*kam.*kap.*kim.*kip.*kma.*kpa.*2.0+cr.*kam.*kap.*kim.*kip.*kma.*kpi.*4.0+cr.*kam.*kap.*kim.*kip.*kmi.*kpa.*4.0+cr.*kam.*kap.*kim.*kip.*kmi.*kpi.*2.0+cr.*kam.*kim.*kip.*kma.*kmi.*kpa.*2.0+cr.*kap.*kim.*kip.*kma.*kmi.*kpa.*2.0+cr.*kam.*kim.*kip.*kma.*kmi.*kpi.*2.0+cr.*kap.*kim.*kip.*kma.*kmi.*kpi.*4.0+kam.*kap.*kim.*kip.*kpa.*kpi.*t2.*2.0+kam.*kim.*kip.*kma.*kpa.*kpi.*t2.*2.0+kap.*kim.*kip.*kma.*kpa.*kpi.*t2.*2.0+kam.*kim.*kip.*kmi.*kpa.*kpi.*t2.*4.0+kap.*kim.*kip.*kmi.*kpa.*kpi.*t2.*2.0+kim.*kip.*kma.*kmi.*kpa.*kpi.*t2.*4.0)+kap.*kim.*kma-kam.*kip.*kmi+kim.*kma.*kmi+cr.*kam.*kip.*kpa-cr.*kap.*kim.*kpi+cr.*kim.*kma.*kpi+cr.*kip.*kmi.*kpa+kip.*kpa.*kpi.*t2).*(-1.0./2.0))./(kim.*kip.*kmi+cr.*kim.*kip.*kpi);
