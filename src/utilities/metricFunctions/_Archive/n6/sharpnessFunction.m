function sharpnessSym = sharpnessFunction(c,k12,k14,k16,k21,k23,k32,k34,k41,k43,k45,k54,k56,k61,k65)
%SHARPNESSFUNCTION
%    SHARPNESSSYM = SHARPNESSFUNCTION(C,K12,K14,K16,K21,K23,K32,K34,K41,K43,K45,K54,K56,K61,K65)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    14-Dec-2020 16:15:23

t2 = c.^2;
t3 = k12.*k14.*k16.*k23.*k45;
t4 = k14.*k16.*k21.*k23.*k45;
t5 = k14.*k16.*k21.*k32.*k45;
t6 = k12.*k14.*k16.*k23.*k65;
t7 = k12.*k14.*k16.*k43.*k45;
t8 = k12.*k16.*k23.*k34.*k45;
t9 = k12.*k16.*k23.*k41.*k45;
t10 = k14.*k16.*k21.*k23.*k65;
t11 = k14.*k16.*k21.*k43.*k45;
t12 = k12.*k16.*k23.*k41.*k54;
t13 = k12.*k16.*k34.*k41.*k45;
t14 = k14.*k16.*k21.*k32.*k65;
t15 = k12.*k14.*k16.*k43.*k65;
t16 = k12.*k14.*k23.*k45.*k56;
t17 = k12.*k16.*k23.*k34.*k65;
t18 = k14.*k16.*k32.*k43.*k45;
t19 = k12.*k14.*k23.*k45.*k61;
t20 = k12.*k16.*k23.*k41.*k65;
t21 = k12.*k16.*k41.*k43.*k45;
t22 = k16.*k21.*k32.*k43.*k45;
t23 = k14.*k16.*k21.*k43.*k65;
t24 = k14.*k21.*k23.*k45.*k56;
t25 = k16.*k23.*k34.*k41.*k45;
t26 = k12.*k14.*k23.*k56.*k61;
t27 = k12.*k16.*k41.*k43.*k54;
t28 = k16.*k21.*k32.*k43.*k54;
t29 = k12.*k16.*k34.*k41.*k65;
t30 = k14.*k21.*k32.*k45.*k56;
t31 = k16.*k32.*k34.*k41.*k45;
t32 = k12.*k14.*k43.*k45.*k56;
t33 = k12.*k16.*k23.*k54.*k65;
t34 = k12.*k23.*k34.*k45.*k56;
t35 = k14.*k16.*k32.*k43.*k65;
t36 = k12.*k14.*k23.*k61.*k65;
t37 = k12.*k14.*k43.*k45.*k61;
t38 = k12.*k23.*k34.*k45.*k61;
t39 = k12.*k16.*k41.*k43.*k65;
t40 = k12.*k23.*k41.*k45.*k56;
t41 = k16.*k21.*k32.*k43.*k65;
t42 = k16.*k32.*k41.*k43.*k45;
t43 = k14.*k21.*k43.*k45.*k56;
t44 = k16.*k21.*k23.*k54.*k65;
t45 = k16.*k23.*k34.*k41.*k65;
t46 = k12.*k14.*k43.*k56.*k61;
t47 = k12.*k23.*k34.*k56.*k61;
t48 = k12.*k23.*k41.*k54.*k56;
t49 = k16.*k32.*k41.*k43.*k54;
t50 = k12.*k34.*k41.*k45.*k56;
t51 = k16.*k21.*k32.*k54.*k65;
t52 = k16.*k32.*k34.*k41.*k65;
t53 = k12.*k16.*k43.*k54.*k65;
t54 = k14.*k32.*k43.*k45.*k56;
t55 = k12.*k14.*k43.*k61.*k65;
t56 = k12.*k23.*k34.*k61.*k65;
t57 = k12.*k23.*k41.*k54.*k65;
t58 = k14.*k32.*k43.*k45.*k61;
t59 = k12.*k23.*k45.*k56.*k61;
t60 = k12.*k41.*k43.*k45.*k56;
t61 = k16.*k32.*k41.*k43.*k65;
t62 = k21.*k32.*k43.*k45.*k56;
t63 = k16.*k21.*k43.*k54.*k65;
t64 = k23.*k34.*k41.*k45.*k56;
t65 = k12.*k23.*k54.*k56.*k61;
t66 = k12.*k41.*k43.*k54.*k56;
t67 = k14.*k32.*k43.*k56.*k61;
t68 = k21.*k32.*k43.*k54.*k56;
t69 = k12.*k34.*k45.*k56.*k61;
t70 = k32.*k34.*k41.*k45.*k56;
t71 = k16.*k32.*k43.*k54.*k65;
t72 = k12.*k23.*k54.*k61.*k65;
t73 = k12.*k41.*k43.*k54.*k65;
t74 = k14.*k32.*k43.*k61.*k65;
t75 = k21.*k32.*k43.*k54.*k65;
t76 = k12.*k43.*k45.*k56.*k61;
t77 = k32.*k41.*k43.*k45.*k56;
t78 = k23.*k34.*k45.*k56.*k61;
t79 = k12.*k43.*k54.*k56.*k61;
t80 = k32.*k41.*k43.*k54.*k56;
t81 = k32.*k34.*k45.*k56.*k61;
t82 = k12.*k43.*k54.*k61.*k65;
t83 = k32.*k41.*k43.*k54.*k65;
t84 = k32.*k43.*k45.*k56.*k61;
t85 = k32.*k43.*k54.*k56.*k61;
t86 = k32.*k43.*k54.*k61.*k65;
t123 = c.*k16.*k21.*k23.*k34.*k45.*2.0;
t124 = c.*k16.*k21.*k32.*k34.*k45.*2.0;
t125 = c.*k16.*k21.*k23.*k34.*k65.*2.0;
t126 = c.*k16.*k21.*k32.*k34.*k65.*2.0;
t127 = c.*k21.*k23.*k34.*k45.*k56.*2.0;
t128 = c.*k21.*k32.*k34.*k45.*k56.*2.0;
t87 = c.*t4;
t88 = c.*t5;
t89 = c.*t8;
t90 = c.*t10;
t91 = c.*t11;
t92 = c.*t13;
t93 = c.*t14;
t94 = c.*t17;
t95 = c.*t22;
t96 = c.*t23;
t97 = c.*t24;
t98 = c.*t25;
t99 = c.*t28;
t100 = c.*t29;
t101 = c.*t30;
t102 = c.*t31;
t103 = c.*t34;
t104 = c.*t38;
t105 = c.*t41;
t106 = c.*t43;
t107 = c.*t44;
t108 = c.*t45;
t109 = c.*t47;
t110 = c.*t50;
t111 = c.*t51;
t112 = c.*t52;
t113 = c.*t56;
t114 = c.*t62;
t115 = c.*t63;
t116 = c.*t64;
t117 = c.*t68;
t118 = c.*t69;
t119 = c.*t70;
t120 = c.*t75;
t121 = c.*t78;
t122 = c.*t81;
t129 = k16.*k21.*k23.*k34.*k45.*t2;
t130 = k16.*k21.*k32.*k34.*k45.*t2;
t131 = k16.*k21.*k23.*k34.*k65.*t2;
t132 = k16.*k21.*k32.*k34.*k65.*t2;
t133 = k21.*k23.*k34.*k45.*k56.*t2;
t134 = k21.*k32.*k34.*k45.*k56.*t2;
t135 = t8+t17+t34;
t136 = t22+t41+t62;
t137 = t28+t47+t68;
t138 = t38+t56+t75;
t143 = t4+t10+t11+t23+t24+t25+t43+t44+t45+t63+t64+t78+t123+t125+t127;
t144 = t5+t13+t14+t29+t30+t31+t50+t51+t52+t69+t70+t81+t124+t126+t128;
t139 = t3+t6+t7+t15+t16+t18+t32+t33+t35+t53+t54+t71+t89+t94+t103;
t140 = t9+t20+t21+t39+t40+t42+t59+t60+t61+t76+t77+t84+t95+t105+t114;
t141 = t12+t26+t27+t46+t48+t49+t65+t66+t67+t79+t80+t85+t99+t109+t117;
t142 = t19+t36+t37+t55+t57+t58+t72+t73+t74+t82+t83+t86+t104+t113+t120;
t147 = t87+t90+t91+t96+t97+t98+t106+t107+t108+t115+t116+t121+t129+t131+t133;
t148 = t88+t92+t93+t100+t101+t102+t110+t111+t112+t118+t119+t122+t130+t132+t134;
t145 = 1.0./t142;
t146 = t145.^2;
t149 = t135.*t145;
t150 = t136.*t145;
t151 = t137.*t145;
t152 = t139.*t145;
t153 = t140.*t145;
t154 = t141.*t145;
t155 = t143.*t145;
t156 = t144.*t145;
t157 = t145.*t147;
t158 = t145.*t148;
t159 = t138.*t139.*t146;
t160 = t138.*t140.*t146;
t161 = t138.*t141.*t146;
t165 = t138.*t146.*t147;
t166 = t138.*t146.*t148;
t169 = t152+t153+t154+t157+t158+1.0;
t162 = -t159;
t163 = -t160;
t164 = -t161;
t167 = -t165;
t168 = -t166;
t170 = 1.0./t169;
t171 = t170.^2;
t172 = t149+t150+t151+t155+t156+t162+t163+t164+t167+t168;
sharpnessSym = t150.*t170+t151.*t170+t156.*t170+t163.*t170+t164.*t170+t168.*t170-t153.*t171.*t172-t154.*t171.*t172-t158.*t171.*t172;
