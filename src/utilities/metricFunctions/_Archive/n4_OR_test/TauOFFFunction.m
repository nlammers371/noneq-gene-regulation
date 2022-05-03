function ETONMean = TauOFFFunction(c,k12,k14,k21,k23,k32,k41,ss3,ss4)
%TAUOFFFUNCTION
%    ETONMEAN = TAUOFFFUNCTION(C,K12,K14,K21,K23,K32,K41,SS3,SS4)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    17-Mar-2021 11:31:20

t2 = c.*k21;
t3 = k12.*k41;
t4 = k32.*k41;
t5 = k32.*t2;
t6 = t3+t4+t5;
t7 = 1.0./t6;
ETONMean = (k14.*ss4.*t7.*(k12+k32+t2)+k23.*ss3.*t7.*(k12+k41+t2))./(k14.*ss4+k23.*ss3);
