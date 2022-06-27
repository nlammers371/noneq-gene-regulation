function out1 = hmFun12(cr,kip,kap,kim,kam,kpi,kmi,kpa,kma)
%HMFUN12
%    OUT1 = HMFUN12(CR,KIP,KAP,KIM,KAM,KPI,KMI,KPA,KMA)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    26-Apr-2021 12:50:33

t2 = kam.^2;
t3 = kap.^2;
t4 = kim.^2;
t5 = kip.^2;
t6 = kma.^2;
t7 = kmi.^2;
t8 = kpa.^2;
t9 = kpi.^2;
out1 = (sqrt(t2.*t3.*t8+t2.*t5.*t8+t3.*t4.*t9+t2.*t7.*t8+t3.*t6.*t9+t4.*t5.*t9+t4.*t6.*t9+t5.*t7.*t8-kam.*kip.*t7.*t8.*2.0-kap.*kim.*t6.*t9.*2.0-kap.*kip.*t2.*t8.*2.0-kap.*kip.*t4.*t9.*2.0-kap.*kma.*t4.*t9.*2.0+kam.*kmi.*t5.*t8.*2.0+kap.*kmi.*t2.*t8.*2.0+kim.*kma.*t3.*t9.*2.0+kip.*kma.*t4.*t9.*2.0-kip.*kmi.*t2.*t8.*2.0-kam.*kap.*kip.*kmi.*t8.*2.0-kap.*kim.*kip.*kma.*t9.*2.0+kam.*kim.*kpa.*kpi.*t3.*2.0+kam.*kim.*kpa.*kpi.*t5.*2.0-kam.*kma.*kpa.*kpi.*t3.*2.0+kam.*kmi.*kpa.*kpi.*t5.*4.0+kim.*kma.*kpa.*kpi.*t3.*4.0-kim.*kmi.*kpa.*kpi.*t5.*2.0-kam.*kap.*kim.*kip.*kpa.*kpi.*4.0-kam.*kap.*kim.*kma.*kpa.*kpi.*2.0+kam.*kap.*kip.*kma.*kpa.*kpi.*2.0+kam.*kap.*kim.*kmi.*kpa.*kpi.*2.0-kam.*kap.*kip.*kmi.*kpa.*kpi.*4.0-kam.*kap.*kma.*kmi.*kpa.*kpi.*2.0+kam.*kim.*kip.*kma.*kpa.*kpi.*2.0-kap.*kim.*kip.*kma.*kpa.*kpi.*4.0-kam.*kim.*kip.*kmi.*kpa.*kpi.*2.0+kap.*kim.*kip.*kmi.*kpa.*kpi.*2.0-kam.*kim.*kma.*kmi.*kpa.*kpi.*2.0+kam.*kip.*kma.*kmi.*kpa.*kpi.*4.0+kap.*kim.*kma.*kmi.*kpa.*kpi.*4.0-kap.*kip.*kma.*kmi.*kpa.*kpi.*2.0-kim.*kip.*kma.*kmi.*kpa.*kpi.*2.0)-kam.*kap.*kpa+kam.*kip.*kpa-kap.*kim.*kpi-kam.*kmi.*kpa-kap.*kma.*kpi+kim.*kip.*kpi+kim.*kma.*kpi+kip.*kmi.*kpa)./(cr.*kap.*kpa.*kpi.*2.0-cr.*kip.*kpa.*kpi.*2.0);