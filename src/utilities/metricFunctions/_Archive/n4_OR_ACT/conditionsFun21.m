function out1 = conditionsFun21(c,k12,k14,k21,k23,k32,k34,k41,k43)
%CONDITIONSFUN21
%    OUT1 = CONDITIONSFUN21(C,K12,K14,K21,K23,K32,K34,K41,K43)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    02-Feb-2021 09:23:07

t2 = c.^2;
t3 = k12.^2;
t4 = k14.^2;
t5 = k21.^2;
t6 = k23.^2;
t7 = k32.^2;
t8 = k34.^2;
t9 = k41.^2;
t10 = k43.^2;
t11 = k12.*k14.*k23;
t12 = k12.*k14.*k43;
t13 = k12.*k23.*k41;
t14 = k14.*k32.*k43;
t15 = k12.*k41.*k43;
t16 = k32.*k41.*k43;
t17 = c.*k12.*k23.*k34;
t18 = c.*k14.*k21.*k43;
t19 = c.*k12.*k34.*k41;
t20 = c.*k21.*k32.*k43;
t21 = t18.*2.0;
t22 = t19.*2.0;
t23 = -t18;
t24 = -t19;
t25 = t11.*t16.*4.0;
t26 = t11.*t20.*2.0;
t27 = c.*k21.*k41.*k43.*t11.*2.0;
t28 = c.*k32.*k34.*k43.*t11.*2.0;
t29 = c.*k21.*k32.*k41.*t12.*4.0;
t30 = c.*k34.*k41.*k43.*t11.*4.0;
t31 = t13.*t20.*2.0;
t32 = c.*k32.*k34.*k41.*t12.*2.0;
t33 = c.*k32.*k34.*k43.*t13.*2.0;
t34 = t12+t15;
t35 = k12.*k23.*k32.*k43.*t4.*2.0;
t36 = k12.*k23.*k32.*k43.*t9.*2.0;
t37 = c.*k12.*k21.*k23.*k43.*t4.*2.0;
t38 = k12.*t4.*t20.*4.0;
t39 = c.*k14.*k23.*k34.*k41.*t3.*2.0;
t40 = c.*k14.*k21.*k32.*k41.*t10.*2.0;
t41 = k43.*t9.*t17.*4.0;
t42 = c.*k12.*k32.*k34.*k43.*t9.*2.0;
t47 = t3.*t4.*t6;
t48 = t3.*t6.*t9;
t49 = t4.*t7.*t10;
t50 = t7.*t9.*t10;
t51 = k14.*k41.*t3.*t6.*2.0;
t52 = k14.*k41.*t7.*t10.*2.0;
t53 = c.*k14.*k34.*t3.*t6.*2.0;
t54 = c.*k14.*k21.*t7.*t10.*2.0;
t55 = c.*k23.*k34.*t3.*t9.*2.0;
t56 = c.*k34.*k41.*t3.*t6.*2.0;
t57 = c.*k21.*k32.*t4.*t10.*2.0;
t58 = c.*k21.*k41.*t7.*t10.*2.0;
t65 = t18+t19;
t66 = k21.*k34.*k43.*t2.*t11.*2.0;
t67 = k21.*k32.*k34.*t2.*t12.*4.0;
t68 = k21.*k34.*k41.*t2.*t12.*2.0;
t69 = k12.*k21.*k23.*k32.*k34.*k43.*t2.*2.0;
t70 = k21.*k34.*k43.*t2.*t13.*4.0;
t71 = k21.*k32.*k34.*t2.*t15.*2.0;
t79 = t2.*t3.*t6.*t8;
t80 = t2.*t4.*t5.*t10;
t81 = t2.*t3.*t8.*t9;
t82 = t2.*t5.*t7.*t10;
t83 = k14.*k32.*t2.*t5.*t10.*2.0;
t84 = k23.*k41.*t2.*t3.*t8.*2.0;
t87 = t13+t16+t20;
t43 = -t27;
t44 = -t29;
t45 = -t30;
t46 = -t32;
t59 = -t37;
t60 = -t38;
t61 = -t39;
t62 = -t40;
t63 = -t41;
t64 = -t42;
t72 = -t55;
t73 = -t57;
t74 = -t66;
t75 = -t67;
t76 = -t70;
t77 = -t71;
t78 = t21+t22;
t85 = -t83;
t86 = -t84;
t88 = t34.*t87.*2.0;
t89 = t11+t14+t17+t87;
t90 = t65+t89;
t91 = t25+t26+t27+t28+t29+t30+t31+t32+t33+t35+t36+t37+t38+t41+t42+t47+t48+t49+t50+t51+t52+t53+t54+t56+t58+t61+t62+t66+t67+t68+t69+t70+t71+t72+t73+t79+t80+t81+t82+t85+t86;
t93 = t25+t26+t28+t31+t33+t35+t36+t39+t40+t43+t44+t45+t46+t47+t48+t49+t50+t51+t52+t53+t54+t55+t56+t57+t58+t59+t60+t63+t64+t68+t69+t74+t75+t76+t77+t79+t80+t81+t82+t83+t84;
t92 = sqrt(t91);
t95 = sqrt(t93);
t94 = -t92;
t96 = -t95;
t97 = t65+t92;
t98 = t78+t92;
t99 = t78+t95;
t100 = t89+t92;
t103 = t90+t95;
t105 = (t95 ~= t90);
t115 = t92+t95;
t101 = (t97 ~= t89);
t102 = t23+t24+t100;
t104 = t23+t24+t89+t94;
t106 = (t103 ~= 0.0);
t108 = t90+t96;
t109 = t15.*t103;
t116 = (t98 ~= t95);
t107 = t15.*t102;
t110 = t15.*t104;
t111 = t15.*t108;
t112 = (t88 ~= t109);
t113 = (t88 ~= t110);
t114 = (t88 ~= t111);
out1 = ((((((((t101 & t105) & t106) & t112) & t113) & t114) & t116) & ((((t100 == t65) | (t92 == t99)) | (t115 == t78)) | (t88 == t107))) | (((((((((((t101 & t105) & t106) & t112) & t113) & t114) & t116) & (t100 ~= t65)) & (t92 ~= t99)) & (t115 ~= t78)) & (t88 ~= t107)) & (t95+t98 ~= 0.0)));
