function sharpnessRightSym = sharpnessRightFunction(c,k12,k14,k16,k21,k23,k32,k34,k41,k43,k45,k54,k56,k61,k65)
%SHARPNESSRIGHTFUNCTION
%    SHARPNESSRIGHTSYM = SHARPNESSRIGHTFUNCTION(C,K12,K14,K16,K21,K23,K32,K34,K41,K43,K45,K54,K56,K61,K65)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    20-Jan-2021 12:14:02

t2 = c.^2;
t3 = k12.*k14.*k16.*k23.*k45;
t4 = k14.*k16.*k21.*k32.*k45;
t5 = k12.*k14.*k16.*k23.*k65;
t6 = k12.*k14.*k16.*k43.*k45;
t7 = k12.*k16.*k23.*k41.*k45;
t8 = k12.*k16.*k23.*k41.*k54;
t9 = k12.*k16.*k34.*k41.*k45;
t10 = k14.*k16.*k21.*k32.*k65;
t11 = k12.*k14.*k16.*k43.*k65;
t12 = k12.*k14.*k23.*k45.*k56;
t13 = k14.*k16.*k32.*k43.*k45;
t14 = k12.*k14.*k23.*k45.*k61;
t15 = k12.*k16.*k23.*k41.*k65;
t16 = k12.*k16.*k41.*k43.*k45;
t17 = k12.*k14.*k23.*k56.*k61;
t18 = k12.*k16.*k41.*k43.*k54;
t19 = k12.*k16.*k34.*k41.*k65;
t20 = k14.*k21.*k32.*k45.*k56;
t21 = k16.*k32.*k34.*k41.*k45;
t22 = k12.*k14.*k43.*k45.*k56;
t23 = k12.*k16.*k23.*k54.*k65;
t24 = k14.*k16.*k32.*k43.*k65;
t25 = k12.*k14.*k23.*k61.*k65;
t26 = k12.*k14.*k43.*k45.*k61;
t27 = k12.*k23.*k34.*k45.*k61;
t28 = k12.*k16.*k41.*k43.*k65;
t29 = k12.*k23.*k41.*k45.*k56;
t30 = k16.*k32.*k41.*k43.*k45;
t31 = k12.*k14.*k43.*k56.*k61;
t32 = k12.*k23.*k41.*k54.*k56;
t33 = k16.*k32.*k41.*k43.*k54;
t34 = k12.*k34.*k41.*k45.*k56;
t35 = k16.*k21.*k32.*k54.*k65;
t36 = k16.*k32.*k34.*k41.*k65;
t37 = k12.*k16.*k43.*k54.*k65;
t38 = k14.*k32.*k43.*k45.*k56;
t39 = k12.*k14.*k43.*k61.*k65;
t40 = k12.*k23.*k34.*k61.*k65;
t41 = k12.*k23.*k41.*k54.*k65;
t42 = k14.*k32.*k43.*k45.*k61;
t43 = k12.*k23.*k45.*k56.*k61;
t44 = k12.*k41.*k43.*k45.*k56;
t45 = k16.*k32.*k41.*k43.*k65;
t46 = k12.*k23.*k54.*k56.*k61;
t47 = k12.*k41.*k43.*k54.*k56;
t48 = k14.*k32.*k43.*k56.*k61;
t49 = k12.*k34.*k45.*k56.*k61;
t50 = k32.*k34.*k41.*k45.*k56;
t51 = k16.*k32.*k43.*k54.*k65;
t52 = k12.*k23.*k54.*k61.*k65;
t53 = k12.*k41.*k43.*k54.*k65;
t54 = k14.*k32.*k43.*k61.*k65;
t55 = k21.*k32.*k43.*k54.*k65;
t56 = k12.*k43.*k45.*k56.*k61;
t57 = k32.*k41.*k43.*k45.*k56;
t58 = k12.*k43.*k54.*k56.*k61;
t59 = k32.*k41.*k43.*k54.*k56;
t60 = k32.*k34.*k45.*k56.*k61;
t61 = k12.*k43.*k54.*k61.*k65;
t62 = k32.*k41.*k43.*k54.*k65;
t63 = k32.*k43.*k45.*k56.*k61;
t64 = k32.*k43.*k54.*k56.*k61;
t65 = k32.*k43.*k54.*k61.*k65;
t66 = c.*k14.*k16.*k21.*k23.*k45;
t68 = c.*k12.*k16.*k23.*k34.*k45;
t69 = c.*k14.*k16.*k21.*k23.*k65;
t70 = c.*k14.*k16.*k21.*k43.*k45;
t73 = c.*k12.*k16.*k23.*k34.*k65;
t74 = c.*k16.*k21.*k32.*k43.*k45;
t75 = c.*k14.*k16.*k21.*k43.*k65;
t76 = c.*k14.*k21.*k23.*k45.*k56;
t77 = c.*k16.*k23.*k34.*k41.*k45;
t78 = c.*k16.*k21.*k32.*k43.*k54;
t82 = c.*k12.*k23.*k34.*k45.*k56;
t84 = c.*k16.*k21.*k32.*k43.*k65;
t85 = c.*k14.*k21.*k43.*k45.*k56;
t86 = c.*k16.*k21.*k23.*k54.*k65;
t87 = c.*k16.*k23.*k34.*k41.*k65;
t88 = c.*k12.*k23.*k34.*k56.*k61;
t93 = c.*k21.*k32.*k43.*k45.*k56;
t94 = c.*k16.*k21.*k43.*k54.*k65;
t95 = c.*k23.*k34.*k41.*k45.*k56;
t96 = c.*k21.*k32.*k43.*k54.*k56;
t100 = c.*k23.*k34.*k45.*k56.*k61;
t102 = c.*k16.*k21.*k32.*k34.*k45.*2.0;
t103 = c.*k16.*k21.*k32.*k34.*k65.*2.0;
t104 = c.*k21.*k32.*k34.*k45.*k56.*2.0;
t67 = c.*t4;
t71 = c.*t9;
t72 = c.*t10;
t79 = c.*t19;
t80 = c.*t20;
t81 = c.*t21;
t83 = c.*t27;
t89 = c.*t34;
t90 = c.*t35;
t91 = c.*t36;
t92 = c.*t40;
t97 = c.*t49;
t98 = c.*t50;
t99 = c.*t55;
t101 = c.*t60;
t105 = k16.*k21.*k23.*k34.*k45.*t2;
t106 = k16.*k21.*k32.*k34.*k45.*t2;
t107 = k16.*k21.*k23.*k34.*k65.*t2;
t108 = k16.*k21.*k32.*k34.*k65.*t2;
t109 = k21.*k23.*k34.*k45.*k56.*t2;
t110 = k21.*k32.*k34.*k45.*k56.*t2;
t111 = t27+t40+t55;
t112 = t3+t5+t6+t11+t12+t13+t22+t23+t24+t37+t38+t51+t68+t73+t82;
t113 = t7+t15+t16+t28+t29+t30+t43+t44+t45+t56+t57+t63+t74+t84+t93;
t114 = t8+t17+t18+t31+t32+t33+t46+t47+t48+t58+t59+t64+t78+t88+t96;
t116 = t4+t9+t10+t19+t20+t21+t34+t35+t36+t49+t50+t60+t102+t103+t104;
t115 = t14+t25+t26+t39+t41+t42+t52+t53+t54+t61+t62+t65+t83+t92+t99;
t119 = t66+t69+t70+t75+t76+t77+t85+t86+t87+t94+t95+t100+t105+t107+t109;
t120 = t67+t71+t72+t79+t80+t81+t89+t90+t91+t97+t98+t101+t106+t108+t110;
t117 = 1.0./t115;
t118 = t117.^2;
t121 = t112.*t117;
t122 = t113.*t117;
t123 = t114.*t117;
t124 = t117.*t119;
t125 = t117.*t120;
t126 = t121+t122+t123+t124+t125+1.0;
t127 = 1.0./t126;
sharpnessRightSym = -t125.*t127.^2.*(t116.*t117+t117.*(k14.*k16.*k21.*k23.*k45+k14.*k16.*k21.*k23.*k65+k14.*k16.*k21.*k43.*k45+k14.*k16.*k21.*k43.*k65+k14.*k21.*k23.*k45.*k56+k16.*k23.*k34.*k41.*k45+k14.*k21.*k43.*k45.*k56+k16.*k21.*k23.*k54.*k65+k16.*k23.*k34.*k41.*k65+k16.*k21.*k43.*k54.*k65+k23.*k34.*k41.*k45.*k56+k23.*k34.*k45.*k56.*k61+c.*k16.*k21.*k23.*k34.*k45.*2.0+c.*k16.*k21.*k23.*k34.*k65.*2.0+c.*k21.*k23.*k34.*k45.*k56.*2.0)+t117.*(k12.*k16.*k23.*k34.*k45+k12.*k16.*k23.*k34.*k65+k12.*k23.*k34.*k45.*k56)+t117.*(k16.*k21.*k32.*k43.*k45+k16.*k21.*k32.*k43.*k65+k21.*k32.*k43.*k45.*k56)+t117.*(k16.*k21.*k32.*k43.*k54+k12.*k23.*k34.*k56.*k61+k21.*k32.*k43.*k54.*k56)-t111.*t112.*t118-t111.*t113.*t118-t111.*t114.*t118-t111.*t118.*t119-t111.*t118.*t120)+t116.*t117.*t127-t111.*t118.*t120.*t127;
