function sharpnessSym = sharpnessFunction(cr,cw,a,ki,ka,kb,ku,wab,wib,wua,wba)
%SHARPNESSFUNCTION
%    SHARPNESSSYM = SHARPNESSFUNCTION(CR,CW,A,KI,KA,KB,KU,WAB,WIB,WUA,WBA)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    24-Jun-2022 12:31:13

t2 = cr+cw;
t3 = a.*ku;
t4 = cw.*kb;
t5 = ka+ki;
t6 = ka.*wab;
t7 = ka.*wba;
t8 = ki.*wib;
t9 = ki.*wua;
t10 = wba.*wib;
t11 = wab.*wua;
t12 = a.^2;
t13 = cr.^2;
t14 = ka.^2;
t15 = kb.^2;
t16 = ki.^2;
t17 = wab.^2;
t18 = wua.^2;
t20 = cr.*kb.*wab;
t19 = t4.*2.0;
t21 = ki.*t6;
t22 = ku.*t6;
t23 = t8.*wba;
t24 = t6.*wua;
t25 = ki+t6;
t26 = ku+t6;
t27 = kb.*t3.*wab;
t29 = t10+wab;
t30 = t11+wba;
t31 = t10+wua;
t32 = ki.*t8;
t33 = kb.*t2;
t34 = ka+t4;
t35 = kb.*t6.*wab;
t37 = t6.*t7;
t45 = t3+t6;
t48 = t6+t8;
t28 = kb.*t23;
t36 = ku.*t19.*wba;
t38 = t32.*wab;
t39 = ku+t34;
t41 = kb.*t24.*wba;
t42 = ka+t19;
t43 = ka+t33;
t44 = ku.*t34;
t46 = t34.*wab;
t47 = ki.*t26;
t49 = t26.*wua;
t52 = ka.*ki.*t29;
t53 = t8+t24;
t54 = ka.*ku.*t30;
t55 = ki.*ku.*t31;
t60 = t6.*t34.*wba;
t61 = t20+t26;
t69 = t9.*t45.*wab;
t71 = t27+t35;
t72 = t18.*t22.*t25;
t78 = cr.*t15.*t48.*wba.*2.0;
t80 = t13.*t15.*t48.*wba;
t81 = t7+t9+t23+t24;
t40 = t28.*wab;
t50 = t39.*wba;
t51 = t44.*wba;
t56 = ku.*t42.*wba;
t57 = t43.*wab;
t58 = t43.*wba;
t59 = ku+t46;
t64 = t6.*t43;
t68 = t8+t49;
t70 = t28+t41;
t74 = t3.*t61;
t75 = t5.*t53;
t82 = t71.*wba.*wua;
t85 = t4.*t81;
t88 = t53.*t78;
t90 = t53.*t80;
t101 = t32+t37+t52+t54+t55;
t62 = ki+t50;
t63 = t6.*t50;
t65 = ki+t58;
t66 = ku+t57;
t67 = t6.*t58;
t73 = t7.*t59;
t83 = ku.*t5.*t68;
t84 = t44+t64;
t87 = t4.*t48.*t70;
t89 = t40+t82;
t91 = t21+t56+t60;
t98 = t75+t85;
t103 = kb.*t101;
t76 = t8.*t62;
t77 = t8.*t65;
t79 = t23.*t66;
t86 = t84.*wab;
t92 = t6.*t91;
t93 = t47+t51+t67;
t94 = t22+t36+t47+t63;
t104 = cr.*t103;
t105 = ku.*t68.*t98;
t107 = t78+t103;
t95 = t93.*wua;
t96 = t21+t73+t76;
t97 = t8.*t94;
t99 = t74+t86;
t110 = a.*t3.*t107.*wua;
t113 = t80+t83+t104;
t100 = t8.*t96;
t102 = t99.*wba.*wua;
t106 = t77+t95;
t109 = t92+t97;
t114 = a.*t3.*t113.*wua;
t108 = t4.*t48.*t106;
t111 = t109.*wua;
t112 = t38+t69+t79+t102;
t115 = t72+t100+t111;
t116 = kb.*t115;
t117 = cr.*t116;
t118 = t88+t116;
t119 = a.*t118;
t120 = t90+t105+t117;
t121 = a.*t120;
t122 = t87+t110+t119;
t123 = t108+t114+t121;
t124 = 1.0./t123;
t125 = t124.^2;
sharpnessSym = ka.*t3.*t124.*(t18.*t71+kb.*t8.*t11)+ka.*t4.*t89.*t124-ka.*t4.*t122.*t125.*(t38+t102+t23.*(t3+t57)+t9.*t26.*wab)+a.*ka.*kb.*t112.*t124-ka.*t3.*t122.*t125.*(t18.*t99+t8.^2+t8.*wua.*(t6.*2.0+t33.*wab+ku.*(a+1.0)))+a.*cr.*ka.*kb.*t89.*t124-a.*cr.*ka.*kb.*t112.*t122.*t125;
