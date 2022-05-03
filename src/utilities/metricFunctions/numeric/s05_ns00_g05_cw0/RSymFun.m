function RSym = RSymFun(cr,ki,ka,kp,km,wap,wip,wma,wpa,wmp,wa,wi)
%RSYMFUN
%    RSYM = RSYMFUN(CR,KI,KA,KP,KM,WAP,WIP,WMA,WPA,WMP,WA,WI)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    18-Oct-2021 13:45:33

t2 = cr.*kp;
t3 = ki.*wip;
t4 = km.*wma;
t5 = ka.*5.0;
t6 = wa.^2;
t7 = wa.^3;
t9 = wap.^2;
t10 = wap.^3;
t12 = wap.^5;
t13 = wi.^2;
t14 = wi.^3;
t16 = wip.^2;
t17 = wip.^3;
t19 = wip.^5;
t20 = wma.^2;
t21 = wma.^3;
t23 = wma.^5;
t24 = wmp.^2;
t25 = wmp.^3;
t27 = wpa.^2;
t28 = wpa.^3;
t30 = wpa.^5;
t35 = ka.*wa.*4.0;
t37 = ki.*wi.*2.0;
t38 = km.*wmp.*2.0;
t8 = t6.^2;
t11 = t9.^2;
t15 = t13.^2;
t18 = t16.^2;
t22 = t20.^2;
t26 = t24.^2;
t29 = t27.^2;
t31 = t2.*2.0;
t32 = t2.*3.0;
t33 = t2.*4.0;
t34 = t2.*5.0;
t36 = t5.*wap;
t39 = t2.*wpa;
t41 = t3.*wip;
t42 = t3.*t16;
t43 = t3.*t17;
t45 = t4.*wma;
t46 = t4.*t20;
t47 = t4.*t21;
t53 = t35.*wap;
t54 = t3.*wi.*2.0;
t55 = t4.*wmp.*2.0;
t56 = t2.*t27;
t57 = t2.*t28;
t59 = t2.*t30;
t61 = ka.*t6.*3.0;
t62 = ka.*t7.*2.0;
t63 = t5.*t9;
t64 = t5.*t10;
t66 = t5.*t12;
t67 = ki.*t13.*3.0;
t68 = ki.*t14.*4.0;
t70 = km.*t24.*3.0;
t71 = km.*t25.*4.0;
t91 = t9.*t35;
t92 = t10.*t35;
t94 = t12.*t35;
t97 = t3.*t13.*3.0;
t100 = t3.*t14.*4.0;
t103 = t4.*t24.*3.0;
t106 = t4.*t25.*4.0;
t40 = ka.*t8;
t44 = t3.*t18;
t48 = t4.*t22;
t49 = t31.*wpa;
t50 = t32.*wpa;
t51 = t33.*wpa;
t52 = t34.*wpa;
t58 = t2.*t29;
t65 = t5.*t11;
t69 = ki.*t15.*5.0;
t72 = km.*t26.*5.0;
t73 = t27.*t31;
t74 = t27.*t32;
t75 = t28.*t31;
t76 = t27.*t33;
t77 = t28.*t32;
t78 = t29.*t31;
t79 = t27.*t34;
t80 = t28.*t33;
t81 = t29.*t32;
t82 = t30.*t31;
t83 = t28.*t34;
t84 = t29.*t33;
t85 = t30.*t32;
t86 = t29.*t34;
t87 = t30.*t33;
t88 = t30.*t34;
t89 = t61.*wap;
t90 = t62.*wap;
t93 = t11.*t35;
t95 = t41.*wi.*2.0;
t96 = t42.*wi.*2.0;
t98 = t43.*wi.*2.0;
t101 = t3.*t15.*5.0;
t102 = t45.*wmp.*2.0;
t104 = t46.*wmp.*2.0;
t105 = t47.*wmp.*2.0;
t108 = t4.*t26.*5.0;
t113 = t9.*t61;
t114 = t9.*t62;
t115 = t10.*t61;
t116 = t10.*t62;
t117 = t11.*t61;
t118 = t11.*t62;
t119 = t12.*t61;
t120 = t12.*t62;
t121 = t13.*t41.*3.0;
t122 = t13.*t42.*3.0;
t123 = t13.*t43.*3.0;
t124 = t14.*t41.*4.0;
t126 = t14.*t42.*4.0;
t127 = t14.*t43.*4.0;
t128 = t15.*t41.*5.0;
t130 = t15.*t42.*5.0;
t131 = t15.*t43.*5.0;
t133 = t24.*t45.*3.0;
t134 = t24.*t46.*3.0;
t135 = t25.*t45.*4.0;
t136 = t24.*t47.*3.0;
t137 = t25.*t46.*4.0;
t139 = t26.*t45.*5.0;
t140 = t25.*t47.*4.0;
t141 = t26.*t46.*5.0;
t143 = t26.*t47.*5.0;
t60 = t40.*wap;
t99 = t44.*wi.*2.0;
t107 = t48.*wmp.*2.0;
t109 = t9.*t40;
t110 = t10.*t40;
t111 = t11.*t40;
t112 = t12.*t40;
t125 = t13.*t44.*3.0;
t129 = t14.*t44.*4.0;
t132 = t15.*t44.*5.0;
t138 = t24.*t48.*3.0;
t142 = t25.*t48.*4.0;
t144 = t26.*t48.*5.0;
RSym = reshape([-t5-t34,t34,0.0,0.0,0.0,0.0,t5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,km,-km-t33-ka.*wap.*5.0,t33,0.0,0.0,0.0,0.0,t36,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t38,-t32-t38-ka.*t9.*5.0,t32,0.0,0.0,0.0,0.0,t63,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t70,-t31-t70-ka.*t10.*5.0,t31,0.0,0.0,0.0,0.0,t64,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t71,-t2-t71-ka.*t11.*5.0,t2,0.0,0.0,0.0,0.0,t65,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t72,-t72-ka.*t12.*5.0,0.0,0.0,0.0,0.0,0.0,t66,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,ki,0.0,0.0,0.0,0.0,0.0,-ki-t35-t39.*5.0,t52,0.0,0.0,0.0,0.0,t35,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t3,0.0,0.0,0.0,0.0,t4,-t3-t4-t39.*4.0-ka.*wa.*wap.*4.0,t51,0.0,0.0,0.0,0.0,t53,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t41,0.0,0.0,0.0,0.0,t55,t39.*-3.0-t41-t55-ka.*t9.*wa.*4.0,t50,0.0,0.0,0.0,0.0,t91,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t42,0.0,0.0,0.0,0.0,t103,t39.*-2.0-t42-t103-ka.*t10.*wa.*4.0,t49,0.0,0.0,0.0,0.0,t92,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t43,0.0,0.0,0.0,0.0,t106,-t39-t43-t106-ka.*t11.*wa.*4.0,t39,0.0,0.0,0.0,0.0,t93,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t44,0.0,0.0,0.0,0.0,t108,-t44-t108-ka.*t12.*wa.*4.0,0.0,0.0,0.0,0.0,0.0,t94,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t37,0.0,0.0,0.0,0.0,0.0,-t37-t56.*5.0-t61,t79,0.0,0.0,0.0,0.0,t61,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t54,0.0,0.0,0.0,0.0,t45,-t45-t54-t56.*4.0-ka.*t6.*wap.*3.0,t76,0.0,0.0,0.0,0.0,t89,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t95,0.0,0.0,0.0,0.0,t102,t56.*-3.0-t95-t102-ka.*t6.*t9.*3.0,t74,0.0,0.0,0.0,0.0,t113,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t96,0.0,0.0,0.0,0.0,t133,t56.*-2.0-t96-t133-ka.*t6.*t10.*3.0,t73,0.0,0.0,0.0,0.0,t115,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t98,0.0,0.0,0.0,0.0,t135,-t56-t98-t135-ka.*t6.*t11.*3.0,t56,0.0,0.0,0.0,0.0,t117,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t99,0.0,0.0,0.0,0.0,t139,-t99-t139-ka.*t6.*t12.*3.0,0.0,0.0,0.0,0.0,0.0,t119,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t67,0.0,0.0,0.0,0.0,0.0,t57.*-5.0-t62-t67,t83,0.0,0.0,0.0,0.0,t62,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t97,0.0,0.0,0.0,0.0,t46,-t46-t57.*4.0-t97-ka.*t7.*wap.*2.0,t80,0.0,0.0,0.0,0.0,t90,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t121,0.0,0.0,0.0,0.0,t104,t57.*-3.0-t104-t121-ka.*t7.*t9.*2.0,t77,0.0,0.0,0.0,0.0,t114,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t122,0.0,0.0,0.0,0.0,t134,t57.*-2.0-t122-t134-ka.*t7.*t10.*2.0,t75,0.0,0.0,0.0,0.0,t116,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t123,0.0,0.0,0.0,0.0,t137,-t57-t123-t137-ka.*t7.*t11.*2.0,t57,0.0,0.0,0.0,0.0,t118,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t125,0.0,0.0,0.0,0.0,t141,-t125-t141-ka.*t7.*t12.*2.0,0.0,0.0,0.0,0.0,0.0,t120,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t68,0.0,0.0,0.0,0.0,0.0,-t40-t58.*5.0-t68,t86,0.0,0.0,0.0,0.0,t40,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t100,0.0,0.0,0.0,0.0,t47,-t47-t58.*4.0-t60-t100,t84,0.0,0.0,0.0,0.0,t60,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t124,0.0,0.0,0.0,0.0,t105,t58.*-3.0-t105-t109-t124,t81,0.0,0.0,0.0,0.0,t109,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t126,0.0,0.0,0.0,0.0,t136,t58.*-2.0-t110-t126-t136,t78,0.0,0.0,0.0,0.0,t110,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t127,0.0,0.0,0.0,0.0,t140,-t58-t111-t127-t140,t58,0.0,0.0,0.0,0.0,t111,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t129,0.0,0.0,0.0,0.0,t143,-t112-t129-t143,0.0,0.0,0.0,0.0,0.0,t112,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t69,0.0,0.0,0.0,0.0,0.0,t59.*-5.0-t69,t88,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t101,0.0,0.0,0.0,0.0,t48,-t48-t59.*4.0-t101,t87,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t128,0.0,0.0,0.0,0.0,t107,t59.*-3.0-t107-t128,t85,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t130,0.0,0.0,0.0,0.0,t138,t59.*-2.0-t130-t138,t82,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t131,0.0,0.0,0.0,0.0,t142,-t59-t131-t142,t59,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t132,0.0,0.0,0.0,0.0,t144,-t132-t144],[36,36]);
