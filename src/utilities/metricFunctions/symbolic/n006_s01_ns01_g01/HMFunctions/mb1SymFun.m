function mb1Sym = mb1SymFun(cr,cw,a,ki,ka,kb,ku,wab,wib,wua,wba)
%MB1SYMFUN
%    MB1SYM = MB1SYMFUN(CR,CW,A,KI,KA,KB,KU,WAB,WIB,WUA,WBA)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    24-Jun-2022 12:12:38

t2 = a.^2;
t3 = cr.^2;
t4 = cw.^2;
t5 = ka.^2;
t6 = ka.^3;
t7 = kb.^2;
t8 = ki.^2;
t9 = ki.^3;
t10 = ku.^2;
t11 = ku.^3;
t12 = wab.^2;
t13 = wib.^2;
t14 = wua.^2;
t15 = cw.*ka.*kb.*ki.*ku.*wab.*wua;
t18 = a.*cr.*ka.*kb.*ki.*ku.*wba.*wib;
t19 = a.*cw.*ka.*kb.*ki.*ku.*wba.*wib;
t21 = cw.*ka.*kb.*ki.*ku.*wba.*wib.*wua;
t23 = a.*cr.*ka.*kb.*ki.*ku.*wab.*wib.*wua;
t16 = cw.*ka.*kb.*t8.*wab.*wib;
t17 = cw.*kb.*ku.*t8.*wib.*wua;
t20 = cw.*kb.*t9.*t13;
t22 = a.*cr.*kb.*t9.*t13;
t24 = a.*t15.*wib;
t25 = a.*cr.*ka.*kb.*t8.*wab.*wib;
t26 = a.*cw.*ka.*kb.*t10.*wba.*wua;
t27 = a.*cr.*kb.*ku.*t8.*wib.*wua;
t29 = cw.*kb.*ki.*t5.*wab.*wba.*wib;
t31 = cw.*kb.*ku.*t5.*wab.*wba.*wua;
t32 = a.*cr.*ka.*kb.*ki.*ku.*t14.*wab;
t33 = cr.*ka.*kb.*ki.*ku.*t2.*wab.*wua;
t34 = a.*cw.*ka.*kb.*ki.*ku.*t14.*wab;
t35 = a.*cr.*kb.*ki.*t5.*wab.*wba.*wib;
t37 = a.*cr.*kb.*ku.*t5.*wab.*wba.*wua;
t39 = a.*cw.*kb.*ki.*t10.*wba.*wib.*wua;
t41 = a.*cw.*kb.*ki.*t10.*t14;
t42 = t18.*wab.*wua;
t44 = cw.*ka.*kb.*t8.*t13.*wba;
t45 = cw.*kb.*ki.*t5.*t12.*wua;
t46 = cw.*kb.*t6.*t12.*wba.*wua;
t47 = a.*cr.*ka.*kb.*t8.*t13.*wba;
t48 = a.*cr.*kb.*ki.*t5.*t12.*wua;
t49 = a.*cr.*kb.*ku.*t8.*t13.*wba;
t50 = cr.*ka.*kb.*t2.*t10.*wba.*wua;
t51 = a.*cw.*kb.*ku.*t8.*t13.*wba;
t52 = a.*cw.*ka.*kb.*t10.*t14.*wab;
t53 = cr.*kb.*ku.*t2.*t8.*wib.*wua;
t54 = a.*cr.*kb.*t6.*t12.*wba.*wua;
t55 = cr.*ka.*kb.*ki.*ku.*t2.*wba.*wib.*wua;
t57 = cr.*kb.*ku.*t2.*t5.*wab.*wba.*wua;
t58 = cr.*kb.*ki.*t2.*t10.*wba.*wib.*wua;
t59 = cr.*kb.*ki.*t2.*t10.*t14;
t60 = a.*cr.*kb.*ku.*t5.*t12.*t14;
t61 = cr.*ka.*kb.*t2.*t10.*t14.*wab;
t62 = a.*cw.*kb.*ku.*t5.*t12.*t14;
t28 = a.*t17;
t30 = t16.*wua;
t36 = t25.*wua;
t38 = a.*t31;
t40 = t29.*wua;
t43 = t24.*wba;
t56 = t35.*wua;
mb1Sym = (t15+t16-t17+t18+t19-t20-t21-t22+t23+t24+t25+t26-t27-t28+t29-t30+t31-t32+t33-t34+t35-t36+t37+t38-t39-t40-t41-t42-t43-t44+t45+t46-t47+t48-t49+t50-t51+t52-t53+t54-t55-t56+t57-t58-t59+t60+t61+t62+sqrt(-(a.*ku.*t9.*t13-ka.*t2.*t11.*t14+ki.*t2.*t11.*t14-a.*ka.*ku.*t8.*t13-a.*ku.*t6.*t12.*t14-a.*t5.*t10.*t14.*wab+a.*t8.*t10.*wib.*wua-t2.*t5.*t10.*t14.*wab+t2.*t8.*t10.*wib.*wua+ka.*ki.*t2.*t10.*t14.*wab-ka.*ki.*t2.*t10.*wib.*wua+a.*ki.*ku.*t5.*t12.*t14+a.*ka.*ki.*t10.*t14.*wab-a.*ka.*ki.*t10.*wib.*wua+a.*ka.*ku.*t8.*wab.*wib.*wua.*2.0-a.*ki.*ku.*t5.*wab.*wib.*wua.*2.0).*(t4.*t7.*t8.*t13.*wba.*4.0+a.*t3.*t7.*t8.*t13.*wba.*4.0-t4.*t5.*t7.*t12.*wba.*wua.*4.0+cr.*cw.*t7.*t8.*t13.*wba.*4.0+a.*cr.*cw.*t7.*t8.*t13.*wba.*4.0-cr.*cw.*t5.*t7.*t12.*wba.*wua.*4.0-a.*t3.*t5.*t7.*t12.*wba.*wua.*4.0-ka.*ki.*t4.*t7.*wab.*wba.*wib.*4.0-ka.*ku.*t4.*t7.*wab.*wba.*wua.*4.0+ki.*ku.*t4.*t7.*wba.*wib.*wua.*4.0-a.*cr.*cw.*t5.*t7.*t12.*wba.*wua.*4.0-cr.*cw.*ka.*ki.*t7.*wab.*wba.*wib.*4.0-a.*ka.*ki.*t3.*t7.*wab.*wba.*wib.*4.0-ka.*ku.*t2.*t3.*t7.*wab.*wba.*wua.*4.0+ki.*ku.*t2.*t3.*t7.*wba.*wib.*wua.*4.0+ka.*ki.*t4.*t7.*wab.*wba.*wib.*wua.*4.0-a.*cr.*cw.*ka.*ki.*t7.*wab.*wba.*wib.*4.0-a.*cr.*cw.*ka.*ku.*t7.*wab.*wba.*wua.*8.0+a.*cr.*cw.*ki.*ku.*t7.*wba.*wib.*wua.*8.0+cr.*cw.*ka.*ki.*t7.*wab.*wba.*wib.*wua.*4.0+a.*ka.*ki.*t3.*t7.*wab.*wba.*wib.*wua.*4.0+a.*cr.*cw.*ka.*ki.*t7.*wab.*wba.*wib.*wua.*4.0)+(t15+t16-t17+t18+t19-t20-t21-t22+t23+t24+t25+t26-t27-t28+t29-t30+t31-t32+t33-t34+t35-t36+t37+t38-t39-t40-t41-t42-t43-t44+t45+t46-t47+t48-t49+t50-t51+t52-t53+t54-t55-t56+t57-t58-t59+t60+t61+t62).^2))./(t4.*t7.*t8.*t13.*wba.*2.0+a.*t3.*t7.*t8.*t13.*wba.*2.0-t4.*t5.*t7.*t12.*wba.*wua.*2.0+cr.*cw.*t7.*t8.*t13.*wba.*2.0+a.*cr.*cw.*t7.*t8.*t13.*wba.*2.0-cr.*cw.*t5.*t7.*t12.*wba.*wua.*2.0-a.*t3.*t5.*t7.*t12.*wba.*wua.*2.0-ka.*ki.*t4.*t7.*wab.*wba.*wib.*2.0-ka.*ku.*t4.*t7.*wab.*wba.*wua.*2.0+ki.*ku.*t4.*t7.*wba.*wib.*wua.*2.0-a.*cr.*cw.*t5.*t7.*t12.*wba.*wua.*2.0-cr.*cw.*ka.*ki.*t7.*wab.*wba.*wib.*2.0-a.*ka.*ki.*t3.*t7.*wab.*wba.*wib.*2.0-ka.*ku.*t2.*t3.*t7.*wab.*wba.*wua.*2.0+ki.*ku.*t2.*t3.*t7.*wba.*wib.*wua.*2.0+ka.*ki.*t4.*t7.*wab.*wba.*wib.*wua.*2.0-a.*cr.*cw.*ka.*ki.*t7.*wab.*wba.*wib.*2.0-a.*cr.*cw.*ka.*ku.*t7.*wab.*wba.*wua.*4.0+a.*cr.*cw.*ki.*ku.*t7.*wba.*wib.*wua.*4.0+cr.*cw.*ka.*ki.*t7.*wab.*wba.*wib.*wua.*2.0+a.*ka.*ki.*t3.*t7.*wab.*wba.*wib.*wua.*2.0+a.*cr.*cw.*ka.*ki.*t7.*wab.*wba.*wib.*wua.*2.0);
