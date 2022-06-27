function productionRateSym = productionRateFunction(cr,cw,a,ki,ka,kb,ku,wab,wib,wua,wba)
%PRODUCTIONRATEFUNCTION
%    PRODUCTIONRATESYM = PRODUCTIONRATEFUNCTION(CR,CW,A,KI,KA,KB,KU,WAB,WIB,WUA,WBA)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    24-Jun-2022 12:20:15

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
t17 = wua.^2;
t19 = cr.*kb.*wab;
t18 = t4.*2.0;
t20 = ki.*t6;
t21 = ku.*t6;
t22 = t8.*wba;
t23 = t6.*wua;
t24 = ki+t6;
t25 = ku+t6;
t26 = t10+wab;
t27 = t11+wba;
t28 = t10+wua;
t29 = ki.*t8;
t30 = kb.*t2;
t31 = ka+t4;
t33 = t6.*t7;
t41 = t6+t8;
t32 = ku.*t18.*wba;
t34 = t29.*wab;
t35 = ku+t31;
t36 = ka+t18;
t37 = ka+t30;
t38 = ku.*t31;
t39 = t31.*wab;
t40 = ki.*t25;
t42 = t25.*wua;
t45 = ka.*ki.*t26;
t46 = t8+t23;
t47 = ka.*ku.*t27;
t48 = ki.*ku.*t28;
t53 = t6.*t31.*wba;
t54 = t19+t25;
t61 = t17.*t21.*t24;
t67 = t13.*t15.*t41.*wba;
t68 = t7+t9+t22+t23;
t43 = t35.*wba;
t44 = t38.*wba;
t49 = ku.*t36.*wba;
t50 = t37.*wab;
t51 = t37.*wba;
t52 = ku+t39;
t57 = t6.*t37;
t60 = t8+t42;
t63 = t3.*t54;
t64 = t5.*t46;
t71 = t4.*t68;
t73 = t46.*t67;
t84 = t29+t33+t45+t47+t48;
t55 = ki+t43;
t56 = t6.*t43;
t58 = ki+t51;
t59 = t6.*t51;
t62 = t7.*t52;
t69 = ku.*t5.*t60;
t70 = t38+t57;
t74 = t20+t49+t53;
t81 = t64+t71;
t86 = cr.*kb.*t84;
t65 = t8.*t55;
t66 = t8.*t58;
t72 = t70.*wab;
t75 = t6.*t74;
t76 = t40+t44+t59;
t77 = t21+t32+t40+t56;
t87 = ku.*t60.*t81;
t92 = t67+t69+t86;
t78 = t76.*wua;
t79 = t20+t62+t65;
t80 = t8.*t77;
t82 = t63+t72;
t93 = a.*t3.*t92.*wua;
t83 = t8.*t79;
t85 = t82.*wba.*wua;
t88 = t66+t78;
t90 = t75+t80;
t89 = t4.*t41.*t88;
t91 = t90.*wua;
t94 = t61+t83+t91;
t95 = cr.*kb.*t94;
t96 = t73+t87+t95;
t97 = a.*t96;
t98 = t89+t93+t97;
t99 = 1.0./t98;
productionRateSym = ka.*t3.*t99.*(t17.*t82+t8.^2+t8.*wua.*(t6.*2.0+t30.*wab+ku.*(a+1.0)))+ka.*t4.*t99.*(t34+t85+t22.*(t3+t50)+t9.*t25.*wab)+a.*cr.*ka.*kb.*t99.*(t34+t85+t22.*(ku+t50)+t9.*wab.*(t3+t6));