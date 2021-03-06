function out1 = conditionsFun12(c,k12,k14,k21,k23,k32,k34,k41,k43)
%CONDITIONSFUN12
%    OUT1 = CONDITIONSFUN12(C,K12,K14,K21,K23,K32,K34,K41,K43)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    02-Feb-2021 09:22:19

t2 = k12.^2;
t3 = k14.^2;
t4 = k21.^2;
t5 = k23.^2;
t6 = k32.^2;
t7 = k34.^2;
t8 = k41.^2;
t9 = k43.^2;
t10 = k14.*k21.*k23;
t11 = k14.*k21.*k32;
t12 = k12.*k23.*k34;
t13 = k12.*k23.*k41;
t14 = k14.*k21.*k43;
t15 = k12.*k34.*k41;
t16 = k12.*k41.*k43;
t17 = k21.*k32.*k43;
t18 = k23.*k34.*k41;
t19 = k32.*k34.*k41;
t20 = k32.*k41.*k43;
t21 = c.*k21.*k23.*k34;
t22 = c.*k21.*k32.*k34;
t23 = t12.*2.0;
t24 = t17.*2.0;
t25 = -t12;
t26 = -t17;
t27 = k12.*k32.*k34.*t10.*2.0;
t28 = t10.*t15.*2.0;
t29 = k12.*k34.*k43.*t10.*2.0;
t30 = t11.*t15.*2.0;
t31 = k12.*k34.*k43.*t11.*4.0;
t32 = t14.*t15.*2.0;
t33 = k21.*k32.*k41.*t12.*4.0;
t35 = t10.*t19.*4.0;
t36 = k32.*k34.*k43.*t10.*4.0;
t37 = k21.*k41.*k43.*t12.*4.0;
t38 = k34.*k41.*k43.*t10.*2.0;
t40 = k34.*k41.*k43.*t11.*2.0;
t42 = k12.*k14.*k21.*k34.*t5.*2.0;
t43 = k21.*t5.*t15.*4.0;
t44 = k14.*k21.*k34.*k41.*t5.*2.0;
t45 = k14.*k23.*k32.*k43.*t4.*2.0;
t46 = k32.*t7.*t13.*2.0;
t47 = k14.*k21.*k34.*k41.*t6.*2.0;
t48 = k34.*t6.*t14.*4.0;
t49 = k21.*k34.*k41.*k43.*t6.*2.0;
t56 = t15.*t17.*-2.0;
t57 = t17.*t18.*-2.0;
t58 = t3.*t4.*t5;
t59 = t3.*t4.*t6;
t60 = t2.*t5.*t7;
t61 = t3.*t4.*t9;
t62 = t2.*t7.*t8;
t63 = t4.*t6.*t9;
t64 = t5.*t7.*t8;
t65 = t6.*t7.*t8;
t66 = k23.*k32.*t3.*t4.*2.0;
t67 = k23.*k43.*t3.*t4.*2.0;
t68 = k12.*k23.*t7.*t8.*2.0;
t69 = k12.*k41.*t5.*t7.*2.0;
t70 = k14.*k32.*t4.*t9.*2.0;
t71 = k14.*k43.*t4.*t6.*2.0;
t72 = k23.*k41.*t2.*t7.*2.0;
t73 = k32.*k43.*t3.*t4.*2.0;
t74 = k12.*k32.*t7.*t8.*2.0;
t75 = k23.*k32.*t7.*t8.*2.0;
t82 = t21+t22;
t88 = t13+t16+t20;
t34 = t17.*t23;
t39 = t15.*t24;
t41 = t18.*t24;
t50 = -t27;
t51 = -t29;
t52 = -t31;
t53 = -t33;
t54 = -t36;
t55 = -t37;
t76 = -t42;
t77 = -t43;
t78 = -t45;
t79 = -t46;
t80 = -t48;
t81 = -t49;
t83 = t23+t24;
t84 = -t69;
t85 = -t70;
t86 = -t71;
t87 = -t72;
t89 = t82.*t88.*2.0;
t90 = t27+t28+t29+t30+t31+t32+t33+t34+t35+t36+t37+t38+t39+t40+t41+t42+t43+t44+t47+t48+t49+t58+t59+t60+t61+t62+t63+t64+t65+t66+t67+t68+t73+t74+t75+t78+t79+t84+t85+t86+t87;
t92 = t28+t30+t32+t34+t35+t38+t40+t44+t45+t46+t47+t50+t51+t52+t53+t54+t55+t56+t57+t58+t59+t60+t61+t62+t63+t64+t65+t66+t67+t68+t69+t70+t71+t72+t73+t74+t75+t76+t77+t80+t81;
t91 = sqrt(t90);
t94 = sqrt(t92);
t93 = -t91;
t95 = -t94;
t96 = t83+t91;
t97 = t83+t94;
t98 = t10+t11+t14+t15+t18+t19+t25+t26+t91;
t99 = t10+t11+t12+t14+t15+t17+t18+t19+t94;
t109 = t91+t94;
t100 = t10+t11+t14+t15+t18+t19+t25+t26+t93;
t101 = t10+t11+t12+t14+t15+t17+t18+t19+t95;
t102 = c.*t17.*t98;
t103 = c.*t17.*t99;
t110 = (t96 ~= t94);
t104 = c.*t17.*t100;
t105 = c.*t17.*t101;
t106 = (t89 ~= t103);
t107 = (t89 ~= t104);
t108 = (t89 ~= t105);
out1 = (((((t106 & t107) & t108) & t110) & (((t91 == t97) | (t89 == t102)) | (t109 == t83))) | ((((((t106 & t107) & t108) & t110) & (t91 ~= t97)) & (t89 ~= t102)) & (t109 ~= t83)));
