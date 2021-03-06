function ssVecSym = steadyStateVecFunction(cr,ki,ka,kb,ku,wab,wib,wua,wba)
%STEADYSTATEVECFUNCTION
%    SSVECSYM = STEADYSTATEVECFUNCTION(CR,KI,KA,KB,KU,WAB,WIB,WUA,WBA)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    23-Jun-2022 21:09:32

t2 = ki.*wib;
t3 = ku.*wua;
t4 = cr.^2;
t5 = kb.^2;
t6 = ki.^2;
t8 = ka.*wab.*wua;
t9 = 1.0./ka;
t10 = 1.0./ku;
t11 = cr.*kb.*ki.*wab;
t12 = cr.*kb.*ku.*wba;
t13 = cr.*kb.*wab.*wua;
t17 = cr.*ka.*kb.*wab.*wba;
t7 = ki.*t3;
t14 = ki.*t8;
t15 = ki.*t2;
t18 = cr.*kb.*t2.*wba;
t21 = t4.*t5.*wab.*wba;
t22 = t2.*t4.*t5.*wba;
t23 = t2+t3+t8+t13;
t16 = cr.*kb.*t7;
t19 = cr.*kb.*t15;
t20 = ka.*t18;
t24 = 1.0./t23;
t25 = t7+t14+t15+t18;
t26 = t11+t12+t17+t21;
t27 = t16+t19+t20+t22;
t28 = t9.*t24.*t25;
t29 = t10.*t24.*t26;
t30 = t9.*t10.*t24.*t27;
t31 = t28+t29+t30+1.0;
t32 = 1.0./t31;
ssVecSym = [t28.*t32,t30.*t32,t29.*t32,t32];
