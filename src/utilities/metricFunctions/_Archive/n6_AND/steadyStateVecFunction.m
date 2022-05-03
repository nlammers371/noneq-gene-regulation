function ssVecSym = steadyStateVecFunction(c,k12,k14,k16,k21,k23,k32,k34,k41,k43,k45,k54,k56,k61,k65)
%STEADYSTATEVECFUNCTION
%    SSVECSYM = STEADYSTATEVECFUNCTION(C,K12,K14,K16,K21,K23,K32,K34,K41,K43,K45,K54,K56,K61,K65)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    20-Jan-2021 12:10:54

t2 = c.^2;
t3 = k12.*k14.*k16.*k23.*k45;
t4 = k12.*k14.*k16.*k23.*k65;
t5 = k12.*k14.*k16.*k43.*k45;
t6 = k12.*k16.*k23.*k41.*k45;
t7 = k12.*k16.*k23.*k41.*k54;
t8 = k12.*k14.*k16.*k43.*k65;
t9 = k12.*k14.*k23.*k45.*k56;
t10 = k14.*k16.*k32.*k43.*k45;
t11 = k12.*k14.*k23.*k45.*k61;
t12 = k12.*k16.*k23.*k41.*k65;
t13 = k12.*k16.*k41.*k43.*k45;
t14 = k12.*k14.*k23.*k56.*k61;
t15 = k12.*k16.*k41.*k43.*k54;
t16 = k12.*k14.*k43.*k45.*k56;
t17 = k12.*k16.*k23.*k54.*k65;
t18 = k14.*k16.*k32.*k43.*k65;
t19 = k12.*k14.*k23.*k61.*k65;
t20 = k12.*k14.*k43.*k45.*k61;
t21 = k12.*k16.*k41.*k43.*k65;
t22 = k12.*k23.*k41.*k45.*k56;
t23 = k16.*k32.*k41.*k43.*k45;
t24 = k12.*k14.*k43.*k56.*k61;
t25 = k12.*k23.*k41.*k54.*k56;
t26 = k16.*k32.*k41.*k43.*k54;
t27 = k12.*k16.*k43.*k54.*k65;
t28 = k14.*k32.*k43.*k45.*k56;
t29 = k12.*k14.*k43.*k61.*k65;
t30 = k12.*k23.*k41.*k54.*k65;
t31 = k14.*k32.*k43.*k45.*k61;
t32 = k12.*k23.*k45.*k56.*k61;
t33 = k12.*k41.*k43.*k45.*k56;
t34 = k16.*k32.*k41.*k43.*k65;
t35 = k12.*k23.*k54.*k56.*k61;
t36 = k12.*k41.*k43.*k54.*k56;
t37 = k14.*k32.*k43.*k56.*k61;
t38 = k16.*k32.*k43.*k54.*k65;
t39 = k12.*k23.*k54.*k61.*k65;
t40 = k12.*k41.*k43.*k54.*k65;
t41 = k14.*k32.*k43.*k61.*k65;
t42 = k12.*k43.*k45.*k56.*k61;
t43 = k32.*k41.*k43.*k45.*k56;
t44 = k12.*k43.*k54.*k56.*k61;
t45 = k32.*k41.*k43.*k54.*k56;
t46 = k12.*k43.*k54.*k61.*k65;
t47 = k32.*k41.*k43.*k54.*k65;
t48 = k32.*k43.*k45.*k56.*k61;
t49 = k32.*k43.*k54.*k56.*k61;
t50 = k32.*k43.*k54.*k61.*k65;
t51 = c.*k14.*k16.*k21.*k23.*k45;
t52 = c.*k14.*k16.*k21.*k32.*k45;
t53 = c.*k12.*k16.*k23.*k34.*k45;
t54 = c.*k14.*k16.*k21.*k23.*k65;
t55 = c.*k14.*k16.*k21.*k43.*k45;
t56 = c.*k12.*k16.*k34.*k41.*k45;
t57 = c.*k14.*k16.*k21.*k32.*k65;
t58 = c.*k12.*k16.*k23.*k34.*k65;
t59 = c.*k16.*k21.*k32.*k43.*k45;
t60 = c.*k14.*k16.*k21.*k43.*k65;
t61 = c.*k14.*k21.*k23.*k45.*k56;
t62 = c.*k16.*k23.*k34.*k41.*k45;
t63 = c.*k16.*k21.*k32.*k43.*k54;
t64 = c.*k12.*k16.*k34.*k41.*k65;
t65 = c.*k14.*k21.*k32.*k45.*k56;
t66 = c.*k16.*k32.*k34.*k41.*k45;
t67 = c.*k12.*k23.*k34.*k45.*k56;
t68 = c.*k12.*k23.*k34.*k45.*k61;
t69 = c.*k16.*k21.*k32.*k43.*k65;
t70 = c.*k14.*k21.*k43.*k45.*k56;
t71 = c.*k16.*k21.*k23.*k54.*k65;
t72 = c.*k16.*k23.*k34.*k41.*k65;
t73 = c.*k12.*k23.*k34.*k56.*k61;
t74 = c.*k12.*k34.*k41.*k45.*k56;
t75 = c.*k16.*k21.*k32.*k54.*k65;
t76 = c.*k16.*k32.*k34.*k41.*k65;
t77 = c.*k12.*k23.*k34.*k61.*k65;
t78 = c.*k21.*k32.*k43.*k45.*k56;
t79 = c.*k16.*k21.*k43.*k54.*k65;
t80 = c.*k23.*k34.*k41.*k45.*k56;
t81 = c.*k21.*k32.*k43.*k54.*k56;
t82 = c.*k12.*k34.*k45.*k56.*k61;
t83 = c.*k32.*k34.*k41.*k45.*k56;
t84 = c.*k21.*k32.*k43.*k54.*k65;
t85 = c.*k23.*k34.*k45.*k56.*k61;
t86 = c.*k32.*k34.*k45.*k56.*k61;
t87 = k16.*k21.*k23.*k34.*k45.*t2;
t88 = k16.*k21.*k32.*k34.*k45.*t2;
t89 = k16.*k21.*k23.*k34.*k65.*t2;
t90 = k16.*k21.*k32.*k34.*k65.*t2;
t91 = k21.*k23.*k34.*k45.*k56.*t2;
t92 = k21.*k32.*k34.*k45.*k56.*t2;
t93 = t3+t4+t5+t8+t9+t10+t16+t17+t18+t27+t28+t38+t53+t58+t67;
t94 = t6+t12+t13+t21+t22+t23+t32+t33+t34+t42+t43+t48+t59+t69+t78;
t95 = t7+t14+t15+t24+t25+t26+t35+t36+t37+t44+t45+t49+t63+t73+t81;
t96 = t11+t19+t20+t29+t30+t31+t39+t40+t41+t46+t47+t50+t68+t77+t84;
t97 = 1.0./t96;
t98 = t51+t54+t55+t60+t61+t62+t70+t71+t72+t79+t80+t85+t87+t89+t91;
t99 = t52+t56+t57+t64+t65+t66+t74+t75+t76+t82+t83+t86+t88+t90+t92;
t100 = t93.*t97;
t101 = t94.*t97;
t102 = t95.*t97;
t103 = t97.*t98;
t104 = t97.*t99;
t105 = t100+t101+t102+t103+t104+1.0;
t106 = 1.0./t105;
ssVecSym = [t100.*t106,t103.*t106,t104.*t106,t101.*t106,t102.*t106,t106];
