function out1 = TauCycleFunction(c,k12,k14,k21,k23,k32,k34,k41,k43,ss1,ss2,ss3,ss4)
%TAUCYCLEFUNCTION
%    OUT1 = TAUCYCLEFUNCTION(C,K12,K14,K21,K23,K32,K34,K41,K43,SS1,SS2,SS3,SS4)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    17-Mar-2021 11:31:21

t2 = c.*k21;
t3 = c.*k34;
t4 = k14.*k23;
t5 = k12.*k41;
t6 = k14.*k43;
t7 = k32.*k41;
t8 = k32.*t2;
t9 = k23.*t3;
t10 = t4+t6+t9;
t11 = t5+t7+t8;
t12 = 1.0./t10;
t13 = 1.0./t11;
out1 = (k14.*ss4.*t13.*(k12+k32+t2)+k23.*ss3.*t13.*(k12+k41+t2))./(k14.*ss4+k23.*ss3)+(k32.*ss2.*t12.*(k14+k43+t3)+k41.*ss1.*t12.*(k23+k43+t3))./(k32.*ss2+k41.*ss1);
