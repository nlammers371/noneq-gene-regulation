function RSym = RSymFun(cr,ki,ka,kp,km,wap,wip,wma,wpa,wmp,wa,wi)
%RSYMFUN
%    RSYM = RSYMFUN(CR,KI,KA,KP,KM,WAP,WIP,WMA,WPA,WMP,WA,WI)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    18-Oct-2021 13:41:00

t2 = cr.*kp;
t3 = ki.*wip;
t4 = km.*wma;
t5 = ka.*3.0;
t6 = wa.^2;
t7 = wap.^2;
t8 = wap.^3;
t10 = wap.^5;
t11 = wi.^2;
t12 = wip.^2;
t13 = wip.^3;
t15 = wip.^5;
t16 = wma.^2;
t17 = wma.^3;
t18 = wmp.^2;
t19 = wmp.^3;
t21 = wpa.^2;
t22 = wpa.^3;
t27 = ka.*wa.*2.0;
t29 = ki.*wi.*2.0;
t30 = km.*wmp.*2.0;
t9 = t7.^2;
t14 = t12.^2;
t20 = t18.^2;
t23 = t2.*2.0;
t24 = t2.*3.0;
t25 = t2.*4.0;
t26 = t2.*5.0;
t28 = t5.*wap;
t31 = t2.*wpa;
t32 = ka.*t6;
t33 = t3.*wip;
t34 = t3.*t12;
t35 = t3.*t13;
t37 = t4.*wma;
t38 = t4.*t16;
t43 = t27.*wap;
t44 = t3.*wi.*2.0;
t45 = t4.*wmp.*2.0;
t46 = t2.*t21;
t47 = t2.*t22;
t49 = t5.*t7;
t50 = t5.*t8;
t52 = t5.*t10;
t53 = ki.*t11.*3.0;
t54 = km.*t18.*3.0;
t55 = km.*t19.*4.0;
t65 = t7.*t27;
t66 = t8.*t27;
t68 = t10.*t27;
t71 = t3.*t11.*3.0;
t75 = t4.*t18.*3.0;
t77 = t4.*t19.*4.0;
t36 = t3.*t14;
t39 = t23.*wpa;
t40 = t24.*wpa;
t41 = t25.*wpa;
t42 = t26.*wpa;
t48 = t32.*wap;
t51 = t5.*t9;
t56 = km.*t20.*5.0;
t57 = t21.*t23;
t58 = t21.*t24;
t59 = t22.*t23;
t60 = t21.*t25;
t61 = t22.*t24;
t62 = t21.*t26;
t63 = t22.*t25;
t64 = t22.*t26;
t67 = t9.*t27;
t69 = t33.*wi.*2.0;
t70 = t34.*wi.*2.0;
t72 = t35.*wi.*2.0;
t74 = t37.*wmp.*2.0;
t76 = t38.*wmp.*2.0;
t78 = t4.*t20.*5.0;
t79 = t7.*t32;
t80 = t8.*t32;
t81 = t9.*t32;
t82 = t10.*t32;
t83 = t11.*t33.*3.0;
t84 = t11.*t34.*3.0;
t85 = t11.*t35.*3.0;
t87 = t18.*t37.*3.0;
t88 = t18.*t38.*3.0;
t89 = t19.*t37.*4.0;
t90 = t19.*t38.*4.0;
t91 = t20.*t37.*5.0;
t92 = t20.*t38.*5.0;
t73 = t36.*wi.*2.0;
t86 = t11.*t36.*3.0;
RSym = reshape([-t5-t26,t26,0.0,0.0,0.0,0.0,t5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,km,-km-t25-ka.*wap.*3.0,t25,0.0,0.0,0.0,0.0,t28,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t30,-t24-t30-ka.*t7.*3.0,t24,0.0,0.0,0.0,0.0,t49,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t54,-t23-t54-ka.*t8.*3.0,t23,0.0,0.0,0.0,0.0,t50,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t55,-t2-t55-ka.*t9.*3.0,t2,0.0,0.0,0.0,0.0,t51,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t56,-t56-ka.*t10.*3.0,0.0,0.0,0.0,0.0,0.0,t52,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,ki,0.0,0.0,0.0,0.0,0.0,-ki-t27-t31.*5.0,t42,0.0,0.0,0.0,0.0,t27,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t3,0.0,0.0,0.0,0.0,t4,-t3-t4-t31.*4.0-ka.*wa.*wap.*2.0,t41,0.0,0.0,0.0,0.0,t43,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t33,0.0,0.0,0.0,0.0,t45,t31.*-3.0-t33-t45-ka.*t7.*wa.*2.0,t40,0.0,0.0,0.0,0.0,t65,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t34,0.0,0.0,0.0,0.0,t75,t31.*-2.0-t34-t75-ka.*t8.*wa.*2.0,t39,0.0,0.0,0.0,0.0,t66,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t35,0.0,0.0,0.0,0.0,t77,-t31-t35-t77-ka.*t9.*wa.*2.0,t31,0.0,0.0,0.0,0.0,t67,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t36,0.0,0.0,0.0,0.0,t78,-t36-t78-ka.*t10.*wa.*2.0,0.0,0.0,0.0,0.0,0.0,t68,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t29,0.0,0.0,0.0,0.0,0.0,-t29-t32-t46.*5.0,t62,0.0,0.0,0.0,0.0,t32,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t44,0.0,0.0,0.0,0.0,t37,-t37-t44-t46.*4.0-t48,t60,0.0,0.0,0.0,0.0,t48,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t69,0.0,0.0,0.0,0.0,t74,t46.*-3.0-t69-t74-t79,t58,0.0,0.0,0.0,0.0,t79,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t70,0.0,0.0,0.0,0.0,t87,t46.*-2.0-t70-t80-t87,t57,0.0,0.0,0.0,0.0,t80,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t72,0.0,0.0,0.0,0.0,t89,-t46-t72-t81-t89,t46,0.0,0.0,0.0,0.0,t81,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t73,0.0,0.0,0.0,0.0,t91,-t73-t82-t91,0.0,0.0,0.0,0.0,0.0,t82,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t53,0.0,0.0,0.0,0.0,0.0,t47.*-5.0-t53,t64,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t71,0.0,0.0,0.0,0.0,t38,-t38-t47.*4.0-t71,t63,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t83,0.0,0.0,0.0,0.0,t76,t47.*-3.0-t76-t83,t61,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t84,0.0,0.0,0.0,0.0,t88,t47.*-2.0-t84-t88,t59,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t85,0.0,0.0,0.0,0.0,t90,-t47-t85-t90,t47,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t86,0.0,0.0,0.0,0.0,t92,-t86-t92],[24,24]);
