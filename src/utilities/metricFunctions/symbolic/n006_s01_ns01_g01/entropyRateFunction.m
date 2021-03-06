function entropyRateSym = entropyRateFunction(cr,cw,a,ki,ka,kb,ku,wab,wib,wua,wba)
%ENTROPYRATEFUNCTION
%    ENTROPYRATESYM = ENTROPYRATEFUNCTION(CR,CW,A,KI,KA,KB,KU,WAB,WIB,WUA,WBA)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    27-Jun-2022 15:45:06

t2 = a.*cr;
t3 = cr+cw;
t4 = a.*ku;
t5 = cr.*kb;
t6 = cw.*kb;
t7 = cw.*ku;
t8 = ka+ki;
t9 = ka.*wab;
t10 = ka.*wba;
t11 = ki.*wib;
t12 = ki.*wua;
t13 = wba.*wib;
t14 = wab.*wua;
t15 = a+1.0;
t16 = a.^2;
t17 = cr.^2;
t18 = ka.^2;
t19 = kb.^2;
t20 = ki.^2;
t21 = wib.^2;
t22 = wua.^2;
t32 = 1.0./a;
t33 = 1.0./cr;
t34 = 1.0./cw;
t35 = 1.0./ka;
t36 = 1.0./kb;
t37 = 1.0./ki;
t38 = 1.0./ku;
t39 = 1.0./wab;
t40 = 1.0./wba;
t41 = 1.0./wib;
t42 = 1.0./wua;
t23 = t6.*2.0;
t24 = t9.*2.0;
t25 = ku.*t2;
t26 = ku.*t6;
t27 = t5.*wab;
t28 = ki.*t9;
t29 = ku.*t9;
t30 = t11.*wba;
t31 = t9.*wua;
t43 = ki+t9;
t44 = ku+t9;
t45 = t13+wab;
t46 = t14+wba;
t47 = t13+wua;
t48 = ku.*t15;
t49 = ki.*t11;
t50 = cw+t2;
t51 = kb.*t3;
t52 = ka+t5;
t53 = ka+t6;
t55 = t9.*t10;
t58 = t3.*t9;
t62 = t11.^2;
t66 = t4+t9;
t69 = t9+t11;
t100 = t9./t11;
t54 = ku.*t23.*wba;
t56 = t49.*wab;
t57 = ku+t53;
t59 = t51.*wab;
t60 = t51.*wba;
t61 = ka+t23;
t63 = kb.*t50;
t64 = ka+t51;
t65 = ku.*t53;
t67 = t53.*wab;
t68 = ki.*t44;
t70 = t44.*wua;
t71 = t4.*t52;
t73 = t4+t44;
t76 = ka.*ki.*t45;
t77 = t11+t31;
t78 = ka.*ku.*t46;
t79 = ki.*ku.*t47;
t80 = ku+t4+t24;
t86 = t9.*t53.*wba;
t87 = t27+t44;
t88 = t12.*t44.*wab;
t94 = t4.*t22.*t44;
t99 = t12.*t66.*wab;
t101 = 1.0./t100;
t102 = log(t100);
t104 = ku.*t22.*t66;
t105 = t22.*t29.*t43;
t113 = t7+t25+t58;
t116 = t5.^2.*t69.*wba;
t117 = t10+t12+t30+t31;
t118 = t22.*t44.*t66;
t72 = ki+t60;
t74 = t57.*wba;
t75 = t65.*wba;
t81 = ku.*t61.*wba;
t82 = ka+t63;
t83 = t64.*wab;
t84 = t64.*wba;
t85 = ku+t67;
t89 = ki.*t73;
t92 = t9.*t64;
t93 = ki.*t80;
t98 = t11+t70;
t103 = log(t101);
t108 = t4.*t87;
t109 = t8.*t77;
t120 = t24+t48+t59;
t122 = kb.*t113.*wba;
t125 = t6.*t117;
t129 = t77.*t116;
t144 = t49+t55+t76+t78+t79;
t90 = ki+t74;
t91 = t9.*t74;
t95 = ki+t84;
t96 = ku+t83;
t97 = t9.*t84;
t106 = ku.*t82.*wba;
t107 = t10.*t85;
t110 = t4+t83;
t112 = t11.*t72.*wib;
t119 = ku.*t8.*t98;
t124 = t65+t92;
t126 = t11.*t120.*wua;
t128 = t26+t71+t92;
t130 = t28+t81+t86;
t135 = t93+t122;
t139 = t109+t125;
t149 = t5.*t144;
t111 = t11.*t90;
t114 = t11.*t95;
t115 = t30.*t96;
t123 = t30.*t110;
t127 = t124.*wab;
t131 = t128.*wba;
t132 = t9.*t130;
t133 = t68+t75+t97;
t134 = t29+t54+t68+t91;
t140 = t135.*wib.*wua;
t143 = t89+t97+t106;
t151 = ku.*t98.*t139;
t162 = t116+t119+t149;
t121 = t114.*wib;
t136 = t133.*wua;
t137 = t28+t107+t111;
t138 = t11.*t134;
t141 = t108+t127;
t145 = t89+t131;
t148 = t143.*wib.*wua;
t156 = t112+t118+t140;
t163 = a.*t4.*t162.*wua;
t142 = t11.*t137;
t146 = t141.*wba.*wua;
t147 = t22.*t141;
t150 = t145.*wib.*wua;
t152 = t114+t136;
t154 = t132+t138;
t158 = t94+t121+t148;
t153 = t6.*t69.*t152;
t155 = t154.*wua;
t157 = t62+t126+t147;
t159 = t104+t121+t150;
t160 = t56+t99+t115+t146;
t161 = t56+t88+t123+t146;
t164 = t105+t142+t155;
t165 = t5.*t164;
t166 = t129+t151+t165;
t167 = a.*t166;
t168 = t153+t163+t167;
t169 = 1.0./t168;
entropyRateSym = t6.*t28.*t102.*t158.*t169+ka.*t6.*t11.*t103.*t161.*t169+kb.*t2.*t28.*t102.*t159.*t169+kb.*t10.*t25.*t157.*t169.*log(t5.*t38.*t42.*wba)+ka.*ki.*t4.*t156.*t169.*log(ka.*t37)+ka.*ki.*t4.*t157.*t169.*log(ki.*t35)+kb.*ki.*t25.*t156.*t169.*log(t5.*t38)+kb.*ki.*t25.*t159.*t169.*log(ku./t5)+t4.*t6.*t10.*t157.*t169.*log((t6.*t42.*wba)./t4)+ki.*t4.*t6.*t156.*t169.*log(t6./t4)+ki.*t4.*t6.*t158.*t169.*log(t4./t6)+ka.*t4.*t6.*t161.*t169.*wua.*log((t4.*t40.*wua)./t6)+ka.*kb.*t2.*t11.*t103.*t160.*t169+ka.*kb.*t25.*t160.*t169.*wua.*log((ku.*t40.*wua)./t5);
