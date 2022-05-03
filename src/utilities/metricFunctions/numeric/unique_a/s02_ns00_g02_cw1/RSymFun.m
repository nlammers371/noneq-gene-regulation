function RSym = RSymFun(cr,cw,a,km,kp,ka,ki,wmp21,wmp12,wa,wi,wpa1,wma1,wap1,wip1,wpa2,wma2,wap2,wip2)
%RSYMFUN
%    RSYM = RSYMFUN(CR,CW,A,KM,KP,KA,KI,WMP21,WMP12,WA,WI,WPA1,WMA1,WAP1,WIP1,WPA2,WMA2,WAP2,WIP2)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    30-Sep-2021 14:25:44

t2 = a.*km;
t3 = cr.*kp;
t4 = cw.*kp;
t5 = ka.*wa;
t6 = ki.*wip1;
t7 = ki.*wip2;
t8 = km.*wma1;
t9 = km.*wma2;
t10 = km.*wmp12;
t11 = km.*wmp21;
t12 = ka.*2.0;
t13 = wma1.^2;
t14 = wma2.^2;
t15 = wpa1.^2;
t16 = wpa2.^2;
t19 = ki.*wi.*2.0;
t33 = -km;
t40 = ka.*wap1.*-2.0;
t41 = ka.*wap2.*-2.0;
t17 = t12.*wap1;
t18 = t12.*wap2;
t20 = t2.*wma1;
t21 = t2.*wma2;
t22 = t2.*wmp12;
t23 = t2.*wmp21;
t24 = t3.*wpa1;
t25 = t3.*wpa2;
t26 = t4.*wpa1;
t27 = t4.*wpa2;
t28 = t5.*wap1;
t29 = t5.*wap2;
t30 = t6.*wip2;
t31 = t8.*wmp12;
t32 = t9.*wmp21;
t36 = -t2;
t37 = -t3;
t38 = -t4;
t42 = -t6;
t43 = -t7;
t44 = -t10;
t45 = -t11;
t46 = t8.*wma1;
t47 = t9.*wma2;
t49 = t6.*wi.*2.0;
t50 = t7.*wi.*2.0;
t51 = t2.*t13;
t52 = t2.*t14;
t53 = t3.*t15;
t54 = t3.*t16;
t55 = t4.*t15;
t56 = t4.*t16;
t68 = t40.*wap2;
t34 = t20.*wmp12;
t35 = t21.*wmp21;
t39 = t28.*wap2;
t48 = t17.*wap2;
t57 = t31.*wma1;
t58 = t32.*wma2;
t59 = t30.*wi.*2.0;
t60 = -t22;
t61 = -t23;
t62 = -t24;
t63 = -t25;
t64 = -t26;
t65 = -t27;
t66 = -t28;
t67 = -t29;
t69 = -t49;
t70 = -t50;
t71 = -t30;
t72 = -t31;
t73 = -t32;
t78 = t13.*t22;
t79 = t14.*t23;
t80 = t15.*t37;
t81 = t16.*t37;
t82 = t15.*t38;
t83 = t16.*t38;
t74 = -t34;
t75 = -t35;
t76 = -t39;
t77 = -t59;
t84 = -t57;
t85 = -t58;
t86 = t13.*t60;
t87 = t14.*t61;
RSym = reshape([t3.*-2.0-t4.*2.0-t12,t3,t4,t3,0.0,0.0,t4,0.0,0.0,t12,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,km,t33+t37+t38+t41,0.0,0.0,t3,0.0,0.0,t4,0.0,0.0,t18,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t2,0.0,t36+t37+t38+t41,0.0,0.0,t3,0.0,0.0,t4,0.0,0.0,t18,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,km,0.0,0.0,t33+t37+t38+t40,t3,t4,0.0,0.0,0.0,0.0,0.0,0.0,t17,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t10,0.0,t11,t44+t45+t68,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t48,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t10,t23,0.0,t44+t61+t68,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t48,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t2,0.0,0.0,0.0,0.0,0.0,t36+t37+t38+t40,t3,t4,0.0,0.0,0.0,0.0,0.0,0.0,t17,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t22,0.0,0.0,0.0,0.0,t11,t45+t60+t68,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t48,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t22,0.0,0.0,0.0,t23,0.0,t60+t61+t68,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t48,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,ki,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-ki-t5+t62+t63+t64+t65,t25,t27,t24,0.0,0.0,t26,0.0,0.0,t5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t7,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t9,-t9+t43+t62+t64+t67,0.0,0.0,t24,0.0,0.0,t26,0.0,0.0,t29,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t7,0.0,0.0,0.0,0.0,0.0,0.0,t21,0.0,-t21+t43+t62+t64+t67,0.0,0.0,t24,0.0,0.0,t26,0.0,0.0,t29,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t6,0.0,0.0,0.0,0.0,0.0,t8,0.0,0.0,-t8+t42+t63+t65+t66,t25,t27,0.0,0.0,0.0,0.0,0.0,0.0,t28,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t30,0.0,0.0,0.0,0.0,0.0,t31,0.0,t32,t71+t72+t73+t76,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t39,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t30,0.0,0.0,0.0,0.0,0.0,t31,t35,0.0,t71+t72+t75+t76,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t39,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t6,0.0,0.0,t20,0.0,0.0,0.0,0.0,0.0,-t20+t42+t63+t65+t66,t25,t27,0.0,0.0,0.0,0.0,0.0,0.0,t28,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t30,0.0,0.0,t34,0.0,0.0,0.0,0.0,t32,t71+t73+t74+t76,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t39,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t30,0.0,0.0,t34,0.0,0.0,0.0,t35,0.0,t71+t74+t75+t76,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t39,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t19,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t19+t80+t81+t82+t83,t54,t56,t53,0.0,0.0,t55,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t50,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t47,-t47+t70+t80+t82,0.0,0.0,t53,0.0,0.0,t55,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t50,0.0,0.0,0.0,0.0,0.0,0.0,t52,0.0,t70+t80+t82+t14.*t36,0.0,0.0,t53,0.0,0.0,t55,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t49,0.0,0.0,0.0,0.0,0.0,t46,0.0,0.0,-t46+t69+t81+t83,t54,t56,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t59,0.0,0.0,0.0,0.0,0.0,t57,0.0,t58,t77+t84+t85,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t59,0.0,0.0,0.0,0.0,0.0,t57,t79,0.0,t77+t84+t87,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t49,0.0,0.0,t51,0.0,0.0,0.0,0.0,0.0,t69+t81+t83+t13.*t36,t54,t56,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t59,0.0,0.0,t78,0.0,0.0,0.0,0.0,t58,t77+t85+t86,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t59,0.0,0.0,t78,0.0,0.0,0.0,t79,0.0,t77+t86+t87],[27,27]);
