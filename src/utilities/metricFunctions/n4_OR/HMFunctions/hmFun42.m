function out1 = hmFun42(cr,kip,kap,kim,kam,kpi,kmi,kpa,kma)
%HMFUN42
%    OUT1 = HMFUN42(CR,KIP,KAP,KIM,KAM,KPI,KMI,KPA,KMA)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    26-Apr-2021 13:18:33

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
out1 = ((-sqrt(t5.*t6.*t8+t4.*t7.*t9+t4.*t8.*t9+kip.*kma.*t4.*t9.*2.0+t2.*t4.*t7.*t10+t2.*t5.*t6.*t11+t2.*t4.*t9.*t10+t2.*t5.*t8.*t11+t2.^2.*t5.*t10.*t11+kam.*kap.*kim.*kmi.*t8.*2.0+cr.*kip.*kpa.*t4.*t9.*2.0-cr.*kim.*kpi.*t5.*t8.*2.0+cr.*kma.*kpa.*t4.*t9.*2.0-cr.*kma.*kpi.*t5.*t6.*2.0-cr.*kmi.*kpa.*t4.*t7.*2.0+kim.*kma.*t2.*t5.*t11.*2.0-kip.*kmi.*t2.*t4.*t10.*2.0+kim.*kpa.*t3.*t5.*t11.*2.0+kma.*kpa.*t3.*t5.*t11.*2.0+kam.*kap.*kim.*kip.*kma.*kmi.*2.0+cr.*kam.*kap.*kim.*kpi.*t8.*4.0+cr.*kam.*kap.*kmi.*kpi.*t8.*2.0-cr.*kip.*kma.*kmi.*kpa.*t4.*2.0+kam.*kap.*kip.*kmi.*t2.*t10.*4.0+kam.*kap.*kip.*kpi.*t3.*t10.*2.0+kam.*kap.*kmi.*kpi.*t3.*t10.*2.0-kim.*kma.*kpa.*kpi.*t2.*t5.*2.0+cr.*kam.*kap.*kim.*kip.*kma.*kpa.*2.0+cr.*kam.*kap.*kim.*kip.*kma.*kpi.*4.0+cr.*kam.*kap.*kim.*kip.*kmi.*kpa.*4.0+cr.*kam.*kap.*kim.*kip.*kmi.*kpi.*2.0+cr.*kam.*kap.*kim.*kma.*kmi.*kpa.*2.0+cr.*kam.*kap.*kip.*kma.*kmi.*kpa.*4.0+cr.*kam.*kap.*kim.*kma.*kmi.*kpi.*2.0+cr.*kam.*kap.*kip.*kma.*kmi.*kpi.*2.0+kam.*kap.*kim.*kip.*kpa.*kpi.*t2.*2.0+kam.*kap.*kim.*kma.*kpa.*kpi.*t2.*4.0+kam.*kap.*kip.*kma.*kpa.*kpi.*t2.*2.0+kam.*kap.*kim.*kmi.*kpa.*kpi.*t2.*2.0+kam.*kap.*kip.*kmi.*kpa.*kpi.*t2.*2.0+kam.*kap.*kma.*kmi.*kpa.*kpi.*t2.*4.0)-kap.*kim.*kma+kam.*kip.*kmi+kam.*kma.*kmi-cr.*kam.*kip.*kpa+cr.*kap.*kim.*kpi+cr.*kam.*kmi.*kpa+cr.*kap.*kma.*kpi+kap.*kpa.*kpi.*t2).*(-1.0./2.0))./(kam.*kap.*kma+cr.*kam.*kap.*kpa);
