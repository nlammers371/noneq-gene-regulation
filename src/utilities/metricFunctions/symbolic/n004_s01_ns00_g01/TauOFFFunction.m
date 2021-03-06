function ETONMean = TauOFFFunction(cr,ki,ka,kb,ku,wab,wib,wua,wba)
%TAUOFFFUNCTION
%    ETONMEAN = TAUOFFFUNCTION(CR,KI,KA,KB,KU,WAB,WIB,WUA,WBA)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    23-Jun-2022 21:10:12

t2 = cr.*kb;
t3 = ka.*ku;
t4 = ki.*wib;
t5 = ku.*wua;
t6 = cr.^2;
t7 = ka.^2;
t8 = kb.^2;
t9 = ki.^2;
t11 = ka.*wab.*wua;
t12 = 1.0./ka;
t13 = 1.0./ku;
t10 = ki.*t5;
t14 = ka.*t2.*wab;
t15 = ki.*t2.*wab;
t16 = ku.*t2.*wba;
t17 = t2.*wab.*wua;
t18 = ki.*t11;
t19 = t7.*wab;
t20 = ki.*t4;
t23 = t2.*t4.*wba;
t26 = t2.^2.*wab.*wba;
t21 = t2.*t10;
t22 = t14.*wba;
t24 = t2.*t20;
t25 = ka.*t23;
t27 = t2.*t23;
t28 = t3+t14+t19;
t29 = t4+t5+t11+t17;
t32 = t10+t18+t20+t23;
t30 = 1.0./t28;
t31 = 1.0./t29;
t33 = t15+t16+t22+t26;
t34 = t21+t24+t25+t27;
t35 = t12.*t31.*t32;
t36 = t13.*t31.*t33;
t37 = (t31.*t34)./t3;
t38 = t35+t36+t37+1.0;
t39 = 1.0./t38;
ETONMean = (ki.*t30.*t39.*(ku+t2+ka.*wab)+t4.*t30.*t36.*t39.*(ka+ku+t2))./(ki.*t39+t4.*t36.*t39);
