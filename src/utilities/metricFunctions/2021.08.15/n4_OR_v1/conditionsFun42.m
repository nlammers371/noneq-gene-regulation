function out1 = conditionsFun42(c,k12,k14,k21,k23,k32,k34,k41,k43)
%CONDITIONSFUN42
%    OUT1 = CONDITIONSFUN42(C,K12,K14,K21,K23,K32,K34,K41,K43)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    20-Jan-2021 20:49:08

t2 = c.^2;
t3 = c.^3;
t5 = k12.^2;
t6 = k14.^2;
t7 = k21.^2;
t8 = k23.^2;
t9 = k32.^2;
t10 = k34.^2;
t11 = k41.^2;
t12 = k43.^2;
t13 = k12.*k23.*k41;
t14 = k14.*k32.*k43;
t15 = k12.*k41.*k43;
t16 = k32.*k41.*k43;
t17 = c.*k14.*k21.*k32;
t18 = c.*k12.*k34.*k41;
t19 = c.*k21.*k32.*k43;
t20 = c.*k23.*k34.*k41;
t21 = c.*k32.*k34.*k41;
t4 = t2.^2;
t22 = t14.*2.0;
t23 = t20.*2.0;
t24 = -t14;
t25 = -t20;
t26 = k21.*k32.*k34.*t2;
t28 = t13.*t17.*2.0;
t29 = c.*k14.*k32.*k34.*t13.*4.0;
t31 = t13.*t19.*2.0;
t33 = c.*k21.*k23.*k41.*t14.*4.0;
t34 = c.*k32.*k34.*k43.*t13.*4.0;
t36 = k12.*k14.*k32.*k41.*t12.*2.0;
t37 = t13.*t14.*-2.0;
t38 = t14+t20;
t39 = t16+t21;
t40 = c.*k12.*k21.*k32.*k41.*t12.*2.0;
t41 = c.*k12.*k23.*k34.*k43.*t11.*2.0;
t42 = k41.*t12.*t17.*4.0;
t44 = t14.*t18.*-2.0;
t47 = t5.*t8.*t11;
t48 = t6.*t9.*t12;
t49 = t5.*t11.*t12;
t50 = k23.*k43.*t5.*t11.*2.0;
t52 = c.*k12.*k34.*t8.*t11.*2.0;
t53 = c.*k14.*k21.*t9.*t12.*2.0;
t54 = c.*k23.*k34.*t5.*t11.*2.0;
t55 = c.*k21.*k43.*t6.*t9.*2.0;
t56 = c.*k34.*k43.*t5.*t11.*2.0;
t63 = k21.*k34.*k41.*t2.*t14.*4.0;
t64 = k21.*k23.*k34.*t2.*t16.*2.0;
t69 = k12.*k21.*k32.*k41.*t3.*t10.*2.0;
t70 = k32.*t2.*t10.*t13.*4.0;
t71 = k14.*k21.*k34.*k43.*t2.*t9.*2.0;
t72 = k21.*k23.*k32.*k41.*t3.*t10.*2.0;
t76 = t2.*t6.*t7.*t9;
t77 = t2.*t5.*t10.*t11;
t79 = t2.*t7.*t9.*t12;
t80 = t2.*t8.*t10.*t11;
t81 = k14.*k34.*t3.*t7.*t9.*2.0;
t82 = k12.*k23.*t2.*t10.*t11.*2.0;
t83 = k14.*k43.*t2.*t7.*t9.*2.0;
t84 = k34.*k43.*t3.*t7.*t9.*2.0;
t89 = t13+t15+t19;
t27 = t13.*t22;
t30 = c.*k12.*k21.*k41.*t22;
t32 = t18.*t22;
t35 = t20.*t22;
t43 = -t29;
t45 = -t33;
t46 = -t34;
t51 = -t36;
t57 = -t41;
t58 = -t42;
t59 = k12.*k14.*k41.*t26.*2.0;
t60 = t13.*t26.*2.0;
t61 = k14.*k23.*k41.*t26.*2.0;
t62 = t15.*t26.*4.0;
t65 = t22+t23;
t66 = -t52;
t67 = -t53;
t68 = -t55;
t74 = -t63;
t75 = -t64;
t78 = t4.*t7.*t9.*t10;
t85 = -t70;
t86 = -t71;
t87 = -t72;
t88 = -t82;
t90 = t39.*t89.*2.0;
t91 = t17+t18+t26+t89;
t73 = -t61;
t92 = t38+t91;
t93 = t27+t28+t29+t30+t31+t32+t33+t34+t35+t36+t40+t42+t47+t48+t49+t50+t54+t56+t57+t59+t60+t61+t62+t63+t64+t66+t67+t68+t69+t70+t72+t76+t77+t78+t79+t80+t81+t83+t84+t86+t88;
t94 = sqrt(t93);
t95 = t28+t30+t31+t35+t37+t40+t41+t43+t44+t45+t46+t47+t48+t49+t50+t51+t52+t53+t54+t55+t56+t58+t59+t60+t62+t69+t71+t73+t74+t75+t76+t77+t78+t79+t80+t81+t82+t83+t84+t85+t87;
t96 = -t94;
t97 = sqrt(t95);
t99 = t38+t94;
t100 = t65+t94;
t102 = t91+t94;
t98 = -t97;
t101 = t65+t97;
t103 = (t99 ~= t91);
t104 = t24+t25+t102;
t105 = t92+t97;
t106 = t24+t25+t91+t96;
t107 = (t97 ~= t92);
t117 = t94+t97;
t118 = (t100 ~= t97);
t108 = (t105 ~= 0.0);
t109 = t16.*t104;
t110 = t92+t98;
t111 = t16.*t105;
t112 = t16.*t106;
t113 = t16.*t110;
t114 = (t90 ~= t111);
t115 = (t90 ~= t112);
t116 = (t90 ~= t113);
out1 = ((((((((t103 & t107) & t108) & t114) & t115) & t116) & t118) & ((((t102 == t38) | (t117 == t65)) | (t94 == t101)) | (t90 == t109))) | (((((((((((t103 & t107) & t108) & t114) & t115) & t116) & t118) & (t102 ~= t38)) & (t117 ~= t65)) & (t94 ~= t101)) & (t90 ~= t109)) & (t97+t100 ~= 0.0)));
