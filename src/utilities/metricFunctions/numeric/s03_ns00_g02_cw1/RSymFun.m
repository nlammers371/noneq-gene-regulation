function RSym = RSymFun(cr,cw,a,ki,ka,kp,km,wap,wip,wma,wpa,wmp,wa,wi)
%RSYMFUN
%    RSYM = RSYMFUN(CR,CW,A,KI,KA,KP,KM,WAP,WIP,WMA,WPA,WMP,WA,WI)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    17-Feb-2022 12:32:49

t2 = a.*km;
t3 = cr.*kp;
t4 = cw.*kp;
t5 = ka.*wa;
t6 = ki.*wip;
t7 = km.*wma;
t8 = km.*wmp;
t9 = ka.*2.0;
t10 = wap.^2;
t11 = wap.^3;
t12 = wip.^2;
t13 = wip.^3;
t14 = wma.^2;
t15 = wmp.^2;
t16 = wpa.^2;
t22 = ki.*wi.*2.0;
t35 = ka.*wap.*-2.0;
t17 = t3.*2.0;
t18 = t3.*3.0;
t19 = t4.*2.0;
t20 = t4.*3.0;
t21 = t9.*wap;
t23 = t8.*2.0;
t24 = t2.*wma;
t25 = t2.*wmp;
t26 = t3.*wpa;
t27 = t4.*wpa;
t28 = t5.*wap;
t29 = t7.*wmp;
t31 = -t3;
t33 = -t4;
t36 = -t6;
t37 = t6.*wip;
t38 = t6.*t12;
t39 = t7.*wma;
t40 = t8.*wmp;
t46 = t6.*wi.*2.0;
t48 = t2.*t14;
t49 = t2.*t15;
t50 = t3.*t16;
t51 = t4.*t16;
t52 = t5.*t10;
t53 = t5.*t11;
t54 = t7.*t15;
t57 = t9.*t10;
t58 = t9.*t11;
t80 = ka.*t10.*-2.0;
t81 = ka.*t11.*-2.0;
t30 = t24.*wmp;
t32 = -t17;
t34 = -t19;
t41 = t25.*2.0;
t42 = t17.*wpa;
t43 = t18.*wpa;
t44 = t19.*wpa;
t45 = t20.*wpa;
t47 = t29.*2.0;
t55 = t29.*wma;
t59 = t23.*wmp;
t60 = t40.*3.0;
t61 = -t26;
t62 = t26.*-2.0;
t63 = -t27;
t64 = t27.*-2.0;
t65 = -t28;
t66 = -t46;
t67 = t49.*2.0;
t68 = t49.*3.0;
t69 = t16.*t17;
t70 = t16.*t18;
t71 = t16.*t19;
t72 = t16.*t20;
t73 = t37.*wi.*2.0;
t74 = t38.*wi.*2.0;
t75 = t54.*2.0;
t77 = t54.*3.0;
t78 = t15.*t24;
t79 = t14.*t25;
t82 = t36.*wip;
t83 = t12.*t36;
t84 = t16.*t31;
t85 = t50.*-2.0;
t86 = t16.*t33;
t87 = t51.*-2.0;
t88 = -t52;
t89 = -t53;
t92 = t15.*t39;
t98 = t15.*t48;
t56 = t30.*2.0;
t76 = t47.*wma;
t90 = -t73;
t91 = -t74;
t93 = t78.*2.0;
t94 = t14.*t41;
t95 = t78.*3.0;
t96 = t92.*2.0;
t97 = t92.*3.0;
t99 = t98.*2.0;
t100 = t98.*3.0;
RSym = reshape([-t9-t18-t20,t20,0.0,0.0,t18,0.0,0.0,0.0,0.0,0.0,t9,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t2,-t2+t32+t34+t35,t19,0.0,0.0,t17,0.0,0.0,0.0,0.0,0.0,t21,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t41,t31+t33-t41+t80,t4,0.0,0.0,t3,0.0,0.0,0.0,0.0,0.0,t57,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t68,-t68+t81,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t58,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,km,0.0,0.0,0.0,-km+t32+t34+t35,t19,0.0,t17,0.0,0.0,0.0,0.0,0.0,0.0,t21,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t8,0.0,0.0,t25,-t8-t25+t31+t33+t80,t4,0.0,t3,0.0,0.0,0.0,0.0,0.0,0.0,t57,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t40,0.0,0.0,t67,-t40-t67+t81,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t58,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t23,0.0,0.0,-t23+t31+t33+t80,t4,t3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t57,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t59,0.0,t49,t40.*-2.0-t49+t81,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t58,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t60,0.0,-t60+t81,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t58,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,ki,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-ki-t5-t26.*3.0-t27.*3.0,t45,0.0,0.0,t43,0.0,0.0,0.0,0.0,0.0,t5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t6,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t24,-t24+t36+t62+t64+t65,t44,0.0,0.0,t42,0.0,0.0,0.0,0.0,0.0,t28,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t37,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t56,-t56+t61+t63+t82+t88,t27,0.0,0.0,t26,0.0,0.0,0.0,0.0,0.0,t52,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t38,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t95,t83+t89-t95,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t53,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t6,0.0,0.0,0.0,0.0,0.0,t7,0.0,0.0,0.0,-t7+t36+t62+t64+t65,t44,0.0,t42,0.0,0.0,0.0,0.0,0.0,0.0,t28,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t37,0.0,0.0,0.0,0.0,0.0,t29,0.0,0.0,t30,-t29-t30+t61+t63+t82+t88,t27,0.0,t26,0.0,0.0,0.0,0.0,0.0,0.0,t52,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t38,0.0,0.0,0.0,0.0,0.0,t54,0.0,0.0,t93,-t54+t83+t89-t93,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t53,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t37,0.0,0.0,0.0,0.0,0.0,0.0,t47,0.0,0.0,-t47+t61+t63+t82+t88,t27,t26,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t52,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t38,0.0,0.0,0.0,0.0,0.0,0.0,t75,0.0,t78,-t75-t78+t83+t89,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t53,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t38,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t77,0.0,-t77+t83+t89,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t53,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t22,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t22-t50.*3.0-t51.*3.0,t72,0.0,0.0,t70,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t46,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t48,-t48+t66+t85+t87,t71,0.0,0.0,t69,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t73,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t94,t79.*-2.0+t84+t86+t90,t51,0.0,0.0,t50,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t74,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t100,t91-t100,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t46,0.0,0.0,0.0,0.0,0.0,t39,0.0,0.0,0.0,-t39+t66+t85+t87,t71,0.0,t69,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t73,0.0,0.0,0.0,0.0,0.0,t55,0.0,0.0,t79,-t55-t79+t84+t86+t90,t51,0.0,t50,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t74,0.0,0.0,0.0,0.0,0.0,t92,0.0,0.0,t99,t91-t92-t99,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t73,0.0,0.0,0.0,0.0,0.0,0.0,t76,0.0,0.0,t55.*-2.0+t84+t86+t90,t51,t50,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t74,0.0,0.0,0.0,0.0,0.0,0.0,t96,0.0,t98,t91-t96-t98,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t74,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t97,0.0,t91-t97],[30,30]);
