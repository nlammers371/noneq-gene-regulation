function out1 = hmFun21(cr,kip,kap,kim,kam,kpi,kmi,kpa,kma)
%HMFUN21
%    OUT1 = HMFUN21(CR,KIP,KAP,KIM,KAM,KPI,KMI,KPA,KMA)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    26-Apr-2021 12:54:47

t2 = cr.^2;
t3 = kam.^2;
t4 = kap.^2;
t5 = kim.^2;
t6 = kip.^2;
t7 = kma.^2;
t8 = kmi.^2;
t9 = kpa.^2;
t10 = kpi.^2;
out1 = ((sqrt(t3.*t4.*t7+t4.*t5.*t7+t3.*t6.*t8+t5.*t6.*t8-kam.*kim.*t4.*t7.*2.0-kam.*kim.*t6.*t8.*2.0+t2.*t3.*t8.*t9+t2.*t4.*t7.*t10+t2.*t5.*t7.*t10+t2.*t6.*t8.*t9+kap.*kip.*kma.*kmi.*t3.*2.0+kap.*kip.*kma.*kmi.*t5.*2.0-cr.*kam.*kpa.*t6.*t8.*2.0+cr.*kam.*kpi.*t4.*t7.*2.0+cr.*kap.*kpi.*t5.*t7.*2.0+cr.*kim.*kpa.*t6.*t8.*2.0+cr.*kip.*kpa.*t3.*t8.*2.0-cr.*kim.*kpi.*t4.*t7.*2.0-kam.*kip.*t2.*t8.*t9.*2.0-kap.*kim.*t2.*t7.*t10.*2.0-kam.*kap.*kim.*kip.*kma.*kmi.*4.0-cr.*kam.*kap.*kim.*kpi.*t7.*2.0-cr.*kam.*kim.*kip.*kpa.*t8.*2.0-cr.*kap.*kma.*kmi.*kpa.*t3.*2.0+cr.*kap.*kma.*kmi.*kpi.*t5.*4.0+cr.*kip.*kma.*kmi.*kpa.*t3.*4.0-cr.*kip.*kma.*kmi.*kpi.*t5.*2.0+cr.*kam.*kap.*kim.*kma.*kmi.*kpa.*2.0-cr.*kam.*kap.*kip.*kma.*kmi.*kpa.*2.0-cr.*kam.*kap.*kim.*kma.*kmi.*kpi.*4.0+cr.*kam.*kap.*kip.*kma.*kmi.*kpi.*2.0-cr.*kam.*kim.*kip.*kma.*kmi.*kpa.*4.0+cr.*kap.*kim.*kip.*kma.*kmi.*kpa.*2.0+cr.*kam.*kim.*kip.*kma.*kmi.*kpi.*2.0-cr.*kap.*kim.*kip.*kma.*kmi.*kpi.*2.0-kam.*kap.*kma.*kmi.*kpa.*kpi.*t2.*2.0-kam.*kim.*kma.*kmi.*kpa.*kpi.*t2.*2.0+kam.*kip.*kma.*kmi.*kpa.*kpi.*t2.*4.0+kap.*kim.*kma.*kmi.*kpa.*kpi.*t2.*4.0-kap.*kip.*kma.*kmi.*kpa.*kpi.*t2.*2.0-kim.*kip.*kma.*kmi.*kpa.*kpi.*t2.*2.0)+kam.*kap.*kma-kap.*kim.*kma+kam.*kip.*kmi-kim.*kip.*kmi+cr.*kam.*kmi.*kpa+cr.*kap.*kma.*kpi-cr.*kim.*kma.*kpi-cr.*kip.*kmi.*kpa).*(-1.0./2.0))./(kam.*kma.*kmi-kim.*kma.*kmi);
