function varSum = intrinsicVarianceFunction(cr,ki,ka,kb,ku,wab,wib,wua,wba)
%INTRINSICVARIANCEFUNCTION
%    VARSUM = INTRINSICVARIANCEFUNCTION(CR,KI,KA,KB,KU,WAB,WIB,WUA,WBA)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    23-Jun-2022 21:10:09

t2 = ka.*ku;
t3 = ki.*ku;
t4 = ki.*wab;
t5 = ki.*wib;
t6 = ku.*wba;
t7 = ku.*wua;
t8 = cr.^2;
t9 = ka.^2;
t10 = kb.^2;
t11 = ki.^2;
t12 = ku.^2;
t13 = cr.*kb.*ki;
t19 = ka.*wab.*wba;
t20 = ka.*wab.*wua;
t21 = 1.0./cr;
t22 = 1.0./ka;
t24 = 1.0./kb;
t25 = 1.0./ku;
t26 = cr.*ka.*kb.*wab;
t27 = cr.*ka.*kb.*wba;
t32 = cr.*kb.*wab.*wba;
t33 = cr.*kb.*wab.*wua;
t14 = ka.*t4;
t15 = ka.*t5;
t16 = t3.*wib;
t17 = t2.*wua;
t18 = t3.*wua;
t23 = 1.0./t9;
t28 = cr.*kb.*t4;
t29 = cr.*kb.*t5;
t30 = cr.*kb.*t6;
t31 = cr.*kb.*t7;
t36 = t9.*wab;
t37 = ki.*t5;
t38 = ku.*t7;
t40 = cr.*kb.*t19;
t43 = t5.*t13;
t45 = t8.*t10.*wba;
t48 = t4+t6+t19+t32;
t49 = t5+t7+t20+t33;
t34 = t14.*wua;
t35 = t17.*wab;
t39 = cr.*kb.*t18;
t41 = t29.*wba;
t42 = t31.*wab;
t44 = cr.*kb.*t15.*wba;
t46 = t45.*wab;
t47 = t5.*t45;
t50 = 1.0./t48;
t51 = 1.0./t49;
t55 = t2+t3+t13+t14+t26+t36;
t57 = t2+t15+t16+t26+t29+t36;
t59 = t2+t3+t13+t27+t30+t45;
t62 = t15+t16+t17+t26+t29+t31+t36+t38;
t63 = t2+t3+t13+t14+t30+t36+t40+t45;
t52 = t18+t34+t37+t41;
t53 = t16+t35+t38+t42;
t56 = t28+t30+t40+t46;
t60 = t16+t26+t29+t31+t35+t38;
t61 = t39+t43+t44+t47;
t54 = 1.0./t53;
t58 = 1.0./t56;
t64 = t22.*t51.*t52;
t65 = t25.*t51.*t56;
t66 = (t51.*t61)./t2;
t67 = t64+t65+t66+1.0;
t68 = 1.0./t67;
t69 = t68.^2;
t70 = t22.*t55.*t58.*t68;
t71 = t23.*t51.*t52.*t54.*t60.*t68;
t72 = (t51.*t54.*t56.*t57.*t68)./t2;
t73 = t21.*t23.*t24.*t50.*t51.*t52.*t63.*t68;
t74 = t22.*t58.*t59.*t66.*t68;
t75 = t22.*t54.*t62.*t66.*t68;
varSum = t69.*(t71+t72+t75).*2.0+t65.*t69.*(t71+t72+t75-t22.*t54.*t57).*2.0+t65.*t69.*(t70+t73+t74-t22.*t55.*t58).*2.0+(t51.^2.*1.0./t58.^2.*t69.*(t70+t73+t74).*2.0)./t12;
