function entropyRateSym = entropyRateFunction(cr,cw,b,kip,kap,kim,kam,kpi,kmi,kpa,kma)
%ENTROPYRATEFUNCTION
%    ENTROPYRATESYM = ENTROPYRATEFUNCTION(CR,CW,B,KIP,KAP,KIM,KAM,KPI,KMI,KPA,KMA)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    26-Apr-2021 14:18:11

t2 = b.^2;
t3 = cr.^2;
t4 = kap.^2;
t5 = kip.^2;
t6 = kma.^2;
t7 = kmi.^2;
t8 = 1.0./b;
t9 = 1.0./cr;
t10 = 1.0./cw;
t11 = 1.0./kap;
t12 = 1.0./kip;
t13 = 1.0./kma;
t14 = 1.0./kmi;
t15 = 1.0./kpa;
t16 = 1.0./kpi;
t17 = kam.*kap.*kip.*kma.*kpa;
t18 = kam.*kap.*kip.*kmi.*kpa;
t19 = kam.*kap.*kma.*kmi.*kpa;
t20 = kap.*kim.*kip.*kma.*kpi;
t21 = kap.*kim.*kip.*kmi.*kpi;
t22 = kam.*kip.*kma.*kmi.*kpa;
t23 = kap.*kim.*kma.*kmi.*kpi;
t24 = kim.*kip.*kma.*kmi.*kpi;
t35 = cr.*kap.*kip.*kma.*kpa.*kpi;
t36 = cr.*kap.*kip.*kmi.*kpa.*kpi;
t37 = cw.*kap.*kip.*kma.*kpa.*kpi;
t38 = cw.*kap.*kip.*kmi.*kpa.*kpi;
t39 = cw.*kap.*kma.*kmi.*kpa.*kpi;
t40 = cw.*kip.*kma.*kmi.*kpa.*kpi;
t46 = b.*cr.*kap.*kip.*kma.*kmi.*kpa;
t47 = b.*cr.*kap.*kip.*kma.*kmi.*kpi;
t48 = b.*cw.*kap.*kip.*kma.*kmi.*kpa;
t49 = b.*cw.*kap.*kip.*kma.*kmi.*kpi;
t50 = b.*cr.*kap.*kma.*kmi.*kpa.*kpi;
t52 = b.*cr.*kip.*kma.*kmi.*kpa.*kpi;
t65 = b.*kam.*kap.*kip.*kma.*kmi.*2.0;
t66 = b.*kap.*kim.*kip.*kma.*kmi.*2.0;
t25 = kap.*t12;
t26 = kip.*t11;
t29 = kam.*kma.*kpa.*t4;
t30 = kam.*kmi.*kpa.*t5;
t31 = kim.*kma.*kpi.*t4;
t32 = kim.*kmi.*kpi.*t5;
t33 = b.*t19;
t34 = b.*t24;
t41 = b.*cr.*t17;
t42 = b.*cr.*t18;
t44 = b.*cr.*t20;
t45 = b.*cr.*t21;
t53 = b.*kam.*kap.*kmi.*t6;
t54 = b.*kam.*kip.*kma.*t7;
t55 = b.*kap.*kim.*kmi.*t6;
t56 = b.*kam.*kip.*kpa.*t7;
t57 = b.*kap.*kim.*kpi.*t6;
t58 = b.*kam.*kma.*kpa.*t7;
t59 = b.*kim.*kip.*kma.*t7;
t60 = b.*kim.*kmi.*kpi.*t6;
t61 = cr.*kma.*kpa.*kpi.*t4;
t62 = cw.*kma.*kpa.*kpi.*t4;
t63 = cr.*kmi.*kpa.*kpi.*t5;
t64 = cw.*kmi.*kpa.*kpi.*t5;
t73 = b.*cw.*kap.*kmi.*kpi.*t6;
t74 = b.*cw.*kip.*kma.*kpa.*t7;
t75 = b.*cw.*t35;
t76 = b.*cw.*t36;
t77 = b.*cr.*t39;
t78 = b.*cr.*t40;
t79 = b.*kam.*t4.*t6;
t80 = b.*kim.*t4.*t6;
t81 = b.*kam.*t5.*t7;
t82 = b.*kim.*t5.*t7;
t83 = b.*cr.*kpi.*t4.*t6;
t84 = b.*cw.*kpi.*t4.*t6;
t85 = b.*cr.*kpa.*t5.*t7;
t86 = b.*cw.*kpa.*t5.*t7;
t87 = kam.*kap.*kmi.*t2.*t6;
t88 = kam.*kip.*kma.*t2.*t7;
t89 = kap.*kim.*kmi.*t2.*t6;
t90 = kim.*kip.*kma.*t2.*t7;
t93 = cr.*t2.*t19;
t94 = cr.*t2.*t22;
t95 = cr.*t2.*t23;
t96 = b.*kap.*kip.*kma.*kpa.*kpi.*t3;
t97 = b.*kap.*kip.*kmi.*kpa.*kpi.*t3;
t98 = cr.*t2.*t24;
t99 = kam.*t2.*t6.*t7;
t100 = kim.*t2.*t6.*t7;
t101 = cr.*kam.*kma.*kpa.*t2.*t7;
t102 = cr.*kap.*kmi.*kpi.*t2.*t6;
t103 = b.*kma.*kpa.*kpi.*t3.*t4;
t104 = cr.*kip.*kma.*kpa.*t2.*t7;
t105 = cr.*kim.*kmi.*kpi.*t2.*t6;
t106 = b.*kmi.*kpa.*kpi.*t3.*t5;
t107 = kap.*kma.*kmi.*kpa.*kpi.*t2.*t3;
t108 = kip.*kma.*kmi.*kpa.*kpi.*t2.*t3;
t27 = log(t25);
t28 = log(t26);
t43 = cr.*t33;
t51 = cr.*t34;
t67 = b.*cr.*t29;
t68 = cr.*t56;
t69 = b.*cr.*t30;
t70 = cr.*t57;
t71 = b.*cr.*t31;
t72 = b.*cr.*t32;
t91 = b.*cw.*t61;
t92 = b.*cw.*t63;
t109 = t18+t19+t21+t23+t29+t31+t33+t36+t38+t39+t50+t56+t58+t61+t62;
t110 = t17+t20+t22+t24+t30+t32+t34+t35+t37+t40+t52+t57+t60+t63+t64;
t112 = t47+t49+t53+t54+t65+t73+t79+t81+t83+t84+t87+t88+t99+t102;
t113 = t46+t48+t55+t59+t66+t74+t80+t82+t85+t86+t89+t90+t100+t104;
t111 = 1.0./t110;
t114 = t42+t43+t45+t67+t68+t71+t76+t77+t91+t93+t95+t97+t101+t103+t107;
t115 = t41+t44+t51+t69+t70+t72+t75+t78+t92+t94+t96+t98+t105+t106+t108;
t116 = t109.*t111;
t117 = t10.*t111.*t112;
t118 = t10.*t111.*t113;
t119 = t10.*t111.*t114;
t120 = t10.*t111.*t115;
t121 = t116+t117+t118+t119+t120+1.0;
t122 = 1.0./t121;
entropyRateSym = kap.*t27.*t122+b.*kmi.*t122.*log(b.*kmi.*t10.*t16)+kap.*t27.*t120.*t122+kip.*t28.*t116.*t122+kip.*t28.*t119.*t122+kam.*t118.*t122.*log(kam./kim)+kim.*t117.*t122.*log(kim./kam)+kma.*t119.*t122.*log(kma.*t9.*t15)+kmi.*t120.*t122.*log(kmi.*t9.*t16)+cr.*kpa.*t117.*t122.*log(cr.*kpa.*t13)+cr.*kpi.*t118.*t122.*log(cr.*kpi.*t14)+b.*kma.*t116.*t122.*log(b.*kma.*t10.*t15)+kpa.*t111.*t112.*t122.*log(cw.*kpa.*t8.*t13)+kpi.*t111.*t113.*t122.*log(cw.*kpi.*t8.*t14);
