function varSum = intrinsicVarianceFunction(c,k12,k14,k16,k21,k23,k32,k34,k41,k43,k45,k54,k56,k61,k65)
%INTRINSICVARIANCEFUNCTION
%    VARSUM = INTRINSICVARIANCEFUNCTION(C,K12,K14,K16,K21,K23,K32,K34,K41,K43,K45,K54,K56,K61,K65)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    19-Jan-2021 12:14:05

t2 = c.^2;
t3 = 1.0./c;
t4 = k12.*k14.*k16.*k23;
t5 = k12.*k14.*k16.*k43;
t6 = k12.*k14.*k16.*k45;
t7 = k12.*k16.*k23.*k41;
t8 = k12.*k14.*k23.*k56;
t9 = k12.*k16.*k23.*k54;
t10 = k14.*k16.*k32.*k43;
t11 = k12.*k14.*k16.*k65;
t12 = k14.*k16.*k32.*k45;
t13 = k12.*k14.*k23.*k61;
t14 = k12.*k16.*k41.*k43;
t15 = k12.*k16.*k41.*k45;
t16 = k12.*k16.*k41.*k54;
t17 = k12.*k14.*k43.*k56;
t18 = k12.*k16.*k43.*k54;
t19 = k12.*k14.*k45.*k56;
t20 = k14.*k16.*k32.*k65;
t21 = k12.*k14.*k43.*k61;
t22 = k12.*k23.*k41.*k54;
t23 = k12.*k14.*k45.*k61;
t24 = k12.*k23.*k41.*k56;
t25 = k16.*k32.*k41.*k43;
t26 = k12.*k16.*k41.*k65;
t27 = k16.*k23.*k41.*k54;
t28 = k16.*k32.*k41.*k45;
t29 = k12.*k14.*k56.*k61;
t30 = k16.*k32.*k41.*k54;
t31 = k12.*k23.*k54.*k56;
t32 = k14.*k32.*k43.*k56;
t33 = k16.*k32.*k43.*k54;
t34 = k12.*k16.*k54.*k65;
t35 = k14.*k32.*k45.*k56;
t36 = k12.*k23.*k54.*k61;
t37 = k12.*k41.*k43.*k54;
t38 = k14.*k32.*k43.*k61;
t39 = k12.*k14.*k61.*k65;
t40 = k12.*k23.*k56.*k61;
t41 = k12.*k41.*k43.*k56;
t42 = k14.*k32.*k45.*k61;
t43 = k12.*k41.*k45.*k56;
t44 = k14.*k23.*k56.*k61;
t45 = k16.*k32.*k41.*k65;
t46 = k16.*k41.*k43.*k54;
t47 = k12.*k41.*k54.*k56;
t48 = k14.*k32.*k56.*k61;
t49 = k12.*k43.*k54.*k56;
t50 = k16.*k32.*k54.*k65;
t51 = k12.*k43.*k54.*k61;
t52 = k32.*k41.*k43.*k54;
t53 = k12.*k41.*k54.*k65;
t54 = k12.*k43.*k56.*k61;
t55 = k14.*k32.*k61.*k65;
t56 = k32.*k41.*k43.*k56;
t57 = k12.*k45.*k56.*k61;
t58 = k14.*k43.*k56.*k61;
t59 = k23.*k41.*k54.*k56;
t60 = k32.*k41.*k45.*k56;
t61 = k12.*k54.*k56.*k61;
t62 = k32.*k41.*k54.*k56;
t63 = k32.*k43.*k54.*k56;
t64 = k32.*k43.*k54.*k61;
t65 = k12.*k54.*k61.*k65;
t66 = k32.*k41.*k54.*k65;
t67 = k32.*k43.*k56.*k61;
t68 = k23.*k54.*k56.*k61;
t69 = k32.*k45.*k56.*k61;
t70 = k41.*k43.*k54.*k56;
t71 = k32.*k54.*k56.*k61;
t72 = k32.*k54.*k61.*k65;
t73 = k43.*k54.*k56.*k61;
t74 = c.*k14.*k16.*k21.*k23;
t75 = c.*k14.*k16.*k21.*k32;
t76 = c.*k12.*k16.*k23.*k34;
t77 = c.*k14.*k16.*k21.*k43;
t78 = c.*k14.*k16.*k21.*k45;
t79 = c.*k12.*k16.*k34.*k41;
t80 = c.*k12.*k16.*k34.*k45;
t81 = c.*k14.*k21.*k32.*k45;
t82 = c.*k16.*k21.*k32.*k43;
t83 = c.*k14.*k21.*k23.*k56;
t84 = c.*k16.*k21.*k23.*k54;
t85 = c.*k16.*k21.*k32.*k45;
t86 = c.*k16.*k23.*k34.*k41;
t87 = c.*k14.*k16.*k21.*k65;
t88 = c.*k14.*k21.*k32.*k56;
t89 = c.*k16.*k21.*k32.*k54;
t90 = c.*k16.*k32.*k34.*k41;
t91 = c.*k12.*k23.*k34.*k56;
t92 = c.*k12.*k16.*k34.*k65;
t93 = c.*k16.*k32.*k34.*k45;
t94 = c.*k12.*k23.*k34.*k61;
t95 = c.*k12.*k34.*k41.*k45;
t96 = c.*k14.*k21.*k32.*k65;
t97 = c.*k14.*k21.*k43.*k56;
t98 = c.*k16.*k21.*k32.*k65;
t99 = c.*k16.*k21.*k43.*k54;
t100 = c.*k14.*k21.*k45.*k56;
t101 = c.*k16.*k34.*k41.*k45;
t102 = c.*k12.*k34.*k41.*k56;
t103 = c.*k12.*k34.*k45.*k56;
t104 = c.*k16.*k32.*k34.*k65;
t105 = c.*k21.*k32.*k43.*k54;
t106 = c.*k12.*k34.*k41.*k65;
t107 = c.*k12.*k34.*k45.*k61;
t108 = c.*k21.*k32.*k43.*k56;
t109 = c.*k32.*k34.*k41.*k45;
t110 = c.*k21.*k23.*k54.*k56;
t111 = c.*k21.*k32.*k45.*k56;
t112 = c.*k23.*k34.*k41.*k56;
t113 = c.*k16.*k21.*k54.*k65;
t114 = c.*k16.*k34.*k41.*k65;
t115 = c.*k12.*k34.*k56.*k61;
t116 = c.*k21.*k32.*k54.*k56;
t117 = c.*k32.*k34.*k41.*k56;
t118 = c.*k32.*k34.*k45.*k56;
t119 = c.*k12.*k34.*k61.*k65;
t120 = c.*k21.*k32.*k54.*k65;
t121 = c.*k32.*k34.*k41.*k65;
t122 = c.*k32.*k34.*k45.*k61;
t123 = c.*k21.*k43.*k54.*k56;
t124 = c.*k23.*k34.*k56.*k61;
t125 = c.*k34.*k41.*k45.*k56;
t126 = c.*k32.*k34.*k56.*k61;
t127 = c.*k32.*k34.*k61.*k65;
t128 = c.*k34.*k45.*k56.*k61;
t129 = k45.*t4;
t130 = k21.*t12;
t131 = k65.*t4;
t132 = k45.*t5;
t133 = k45.*t7;
t134 = k54.*t7;
t135 = k34.*t15;
t136 = k21.*t20;
t137 = k65.*t5;
t138 = k45.*t8;
t139 = k45.*t10;
t140 = k45.*t13;
t141 = k65.*t7;
t142 = k45.*t14;
t143 = k61.*t8;
t144 = k54.*t14;
t145 = k34.*t26;
t146 = k21.*t35;
t147 = k34.*t28;
t148 = k45.*t17;
t149 = k65.*t9;
t150 = k65.*t10;
t151 = k65.*t13;
t152 = k45.*t21;
t153 = k65.*t14;
t154 = k45.*t24;
t155 = k45.*t25;
t156 = k61.*t17;
t157 = k56.*t22;
t158 = k54.*t25;
t159 = k34.*t43;
t160 = k21.*t50;
t161 = k34.*t45;
t162 = k65.*t18;
t163 = k45.*t32;
t164 = k65.*t21;
t165 = k65.*t22;
t166 = k45.*t38;
t167 = k45.*t40;
t168 = k45.*t41;
t169 = k65.*t25;
t170 = k61.*t31;
t171 = k56.*t37;
t172 = k61.*t32;
t173 = k34.*t57;
t174 = k34.*t60;
t175 = k65.*t33;
t176 = k65.*t36;
t177 = k65.*t37;
t178 = k65.*t38;
t179 = k45.*t54;
t180 = k45.*t56;
t181 = k61.*t49;
t182 = k56.*t52;
t183 = k34.*t69;
t184 = k65.*t51;
t185 = k65.*t52;
t186 = k45.*t67;
t187 = k61.*t63;
t188 = k65.*t64;
t189 = k45.*t74;
t191 = k45.*t76;
t192 = k65.*t74;
t193 = k45.*t77;
t196 = k34.*t85;
t197 = k65.*t76;
t198 = k45.*t82;
t199 = k65.*t77;
t200 = k45.*t83;
t201 = k45.*t86;
t202 = c.*k21.*t33;
t205 = k34.*t98;
t207 = k45.*t91;
t208 = k45.*t94;
t209 = k65.*t82;
t210 = k45.*t97;
t211 = k65.*t84;
t212 = k65.*t86;
t213 = c.*k34.*t40;
t217 = k34.*t111;
t218 = k65.*t94;
t219 = k45.*t108;
t220 = k65.*t99;
t221 = k45.*t112;
t222 = c.*k21.*t63;
t225 = k65.*t105;
t226 = k45.*t124;
t228 = k16.*k21.*k23.*k34.*t2;
t229 = k16.*k21.*k32.*k34.*t2;
t230 = k16.*k21.*k34.*k45.*t2;
t231 = k21.*k32.*k34.*k45.*t2;
t232 = k21.*k23.*k34.*k56.*t2;
t233 = k16.*k21.*k34.*k65.*t2;
t234 = k21.*k32.*k34.*k56.*t2;
t235 = k21.*k32.*k34.*k65.*t2;
t236 = k21.*k34.*k45.*k56.*t2;
t262 = t6+t11+t12+t15+t16+t19+t20+t23+t26+t28+t29+t30+t34+t35+t39+t42+t43+t45+t47+t48+t50+t53+t55+t57+t60+t61+t62+t65+t66+t69+t71+t72+t78+t85+t87+t89+t98+t100+t111+t113+t116+t120;
t190 = c.*t130;
t194 = c.*t135;
t195 = c.*t136;
t203 = c.*t145;
t204 = c.*t146;
t206 = c.*t147;
t214 = c.*t159;
t215 = c.*t160;
t216 = c.*t161;
t223 = c.*t173;
t224 = c.*t174;
t227 = c.*t183;
t237 = k45.*t228;
t238 = k45.*t229;
t239 = k65.*t228;
t240 = k65.*t229;
t241 = k45.*t232;
t242 = k56.*t231;
t243 = t129+t131+t132+t137+t138+t139+t148+t149+t150+t162+t163+t175+t191+t197+t207;
t244 = t130+t135+t136+t145+t146+t147+t159+t160+t161+t173+t174+t183+t196+t205+t217;
t245 = t133+t141+t142+t153+t154+t155+t167+t168+t169+t179+t180+t186+t198+t209+t219;
t246 = t134+t143+t144+t156+t157+t158+t170+t171+t172+t181+t182+t187+t202+t213+t222;
t247 = t140+t151+t152+t164+t165+t166+t176+t177+t178+t184+t185+t188+t208+t218+t225;
t255 = t4+t5+t7+t9+t10+t13+t14+t18+t21+t22+t25+t33+t36+t37+t38+t51+t52+t64+t74+t75+t76+t77+t79+t82+t84+t86+t89+t90+t94+t99+t105+t228+t229;
t256 = t6+t11+t15+t16+t19+t23+t26+t29+t34+t39+t43+t47+t53+t57+t61+t65+t78+t80+t87+t92+t100+t101+t103+t107+t113+t114+t115+t119+t125+t128+t230+t233+t236;
t263 = t4+t5+t7+t8+t10+t13+t14+t17+t21+t24+t25+t32+t38+t40+t41+t54+t56+t67+t74+t75+t76+t77+t79+t82+t83+t86+t88+t90+t91+t94+t97+t102+t108+t112+t115+t117+t124+t126+t228+t229+t232+t234;
t264 = t6+t11+t12+t15+t16+t19+t20+t23+t26+t28+t29+t30+t34+t35+t39+t42+t43+t45+t47+t48+t50+t53+t55+t57+t60+t61+t62+t65+t66+t69+t71+t72+t78+t80+t87+t92+t93+t100+t103+t104+t107+t113+t115+t118+t119+t122+t126+t127+t230+t233+t236;
t265 = t4+t5+t7+t8+t9+t10+t13+t14+t17+t18+t21+t24+t25+t31+t32+t33+t36+t38+t41+t49+t51+t56+t63+t64+t74+t75+t76+t77+t79+t82+t83+t84+t86+t88+t89+t90+t91+t94+t97+t99+t102+t108+t110+t112+t116+t117+t123+t228+t229+t232+t234;
t266 = t4+t5+t7+t8+t9+t10+t13+t14+t16+t17+t21+t24+t25+t27+t29+t30+t31+t32+t36+t38+t41+t44+t47+t48+t54+t56+t59+t61+t62+t67+t68+t71+t74+t75+t76+t77+t79+t82+t83+t84+t86+t88+t89+t90+t91+t94+t97+t102+t108+t110+t112+t115+t116+t117+t124+t126+t228+t229+t232+t234;
t267 = t6+t11+t12+t15+t16+t19+t20+t23+t26+t28+t29+t30+t34+t35+t39+t42+t43+t45+t47+t48+t50+t53+t55+t57+t60+t61+t62+t65+t66+t69+t71+t72+t75+t78+t79+t85+t87+t88+t89+t90+t92+t96+t100+t102+t104+t106+t111+t113+t115+t116+t117+t119+t120+t121+t126+t127+t229+t233+t234+t235;
t268 = t4+t5+t7+t8+t9+t10+t13+t14+t17+t18+t21+t24+t25+t27+t30+t31+t32+t36+t38+t41+t44+t46+t48+t49+t51+t56+t58+t59+t62+t67+t68+t70+t71+t73+t74+t75+t76+t77+t79+t82+t83+t84+t86+t88+t89+t90+t91+t94+t97+t99+t102+t108+t110+t112+t116+t117+t123+t124+t126+t228+t229+t232+t234;
t269 = t6+t11+t12+t15+t16+t19+t20+t23+t26+t28+t29+t30+t34+t35+t39+t42+t43+t45+t47+t48+t50+t53+t55+t57+t60+t61+t62+t65+t66+t69+t71+t72+t78+t80+t81+t87+t88+t92+t93+t95+t96+t100+t102+t104+t106+t107+t109+t111+t113+t115+t116+t117+t119+t120+t121+t122+t126+t127+t230+t231+t233+t234+t235;
t248 = 1.0./t244;
t249 = 1.0./t246;
t250 = 1.0./t247;
t252 = t189+t192+t193+t199+t200+t201+t210+t211+t212+t220+t221+t226+t237+t239+t241;
t253 = t190+t194+t195+t203+t204+t206+t214+t215+t216+t223+t224+t227+t238+t240+t242;
t251 = t250.^2;
t254 = 1.0./t253;
t257 = t243.*t250;
t258 = t245.*t250;
t259 = t246.*t250;
t260 = t250.*t252;
t261 = t250.*t253;
t270 = t257+t258+t259+t260+t261+1.0;
t271 = 1.0./t270;
t272 = t271.^2;
t273 = t249.*t255.*t271;
t274 = t254.*t269.*t271;
t275 = t254.*t256.*t260.*t271;
t276 = t254.*t258.*t262.*t271;
t277 = t249.*t258.*t263.*t271;
t278 = t3.*t248.*t257.*t264.*t271;
t279 = t249.*t257.*t265.*t271;
t280 = t249.*t261.*t266.*t271;
t281 = t254.*t259.*t267.*t271;
t282 = t249.*t260.*t268.*t271;
varSum = t251.*1.0./t254.^2.*t272.*(t274+t275+t276+t278+t281).*2.0+1.0./t249.^2.*t251.*t272.*(t273+t277+t279+t280+t282).*2.0+t246.*t251.*t253.*t272.*(t274+t275+t276+t278+t281-t254.*t267).*2.0+t246.*t251.*t253.*t272.*(t273+t277+t279+t280+t282-t249.*t266).*2.0;
