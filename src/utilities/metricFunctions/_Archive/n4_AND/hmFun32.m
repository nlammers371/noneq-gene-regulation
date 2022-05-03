function out1 = hmFun32(c,k12,k14,k21,k23,k32,k34,k41,k43)
%HMFUN32
%    OUT1 = HMFUN32(C,K12,K14,K21,K23,K32,K34,K41,K43)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    02-Feb-2021 10:03:56

t2 = c.^2;
t3 = c.^3;
t4 = k12.^2;
t5 = k14.^2;
t6 = k21.^2;
t7 = k23.^2;
t8 = k32.^2;
t9 = k34.^2;
t10 = k41.^2;
t11 = k43.^2;
out1 = ((-sqrt(t4.*t5.*t11+t4.*t7.*t10+t5.*t8.*t11+k12.*k32.*t5.*t11.*2.0+t2.*t5.*t6.*t8+t2.*t4.*t7.*t9+t2.*t5.*t6.*t11+t2.*t7.*t9.*t10+t2.^2.*t6.*t7.*t9-k14.*k23.*k41.*k43.*t4.*2.0+c.*k12.*k21.*t5.*t11.*2.0+c.*k12.*k34.*t7.*t10.*2.0+c.*k21.*k32.*t5.*t11.*2.0-c.*k21.*k43.*t5.*t8.*2.0+c.*k34.*k41.*t4.*t7.*2.0+k12.*k21.*t3.*t7.*t9.*2.0+k12.*k41.*t2.*t7.*t9.*2.0+k21.*k41.*t3.*t7.*t9.*2.0-k32.*k43.*t2.*t5.*t6.*2.0-k12.*k14.*k23.*k32.*k41.*k43.*2.0-c.*k12.*k21.*k32.*k43.*t5.*2.0+c.*k14.*k23.*k34.*k41.*t4.*4.0+c.*k14.*k23.*k34.*k43.*t4.*2.0+k14.*k23.*k32.*k34.*t3.*t6.*2.0+k12.*k21.*k34.*k41.*t2.*t7.*2.0-k14.*k23.*k32.*k43.*t2.*t6.*4.0+k14.*k23.*k34.*k43.*t3.*t6.*2.0-c.*k12.*k14.*k21.*k23.*k32.*k41.*2.0-c.*k12.*k14.*k21.*k23.*k32.*k43.*4.0-c.*k12.*k14.*k21.*k23.*k41.*k43.*2.0+c.*k12.*k14.*k23.*k32.*k34.*k41.*4.0+c.*k12.*k14.*k23.*k32.*k34.*k43.*2.0+c.*k12.*k14.*k23.*k34.*k41.*k43.*2.0-c.*k14.*k21.*k23.*k32.*k41.*k43.*4.0+c.*k14.*k23.*k32.*k34.*k41.*k43.*2.0+k12.*k14.*k21.*k23.*k32.*k34.*t2.*2.0+k12.*k14.*k21.*k23.*k34.*k41.*t2.*4.0+k12.*k14.*k21.*k23.*k34.*k43.*t2.*4.0+k14.*k21.*k23.*k32.*k34.*k41.*t2.*2.0+k14.*k21.*k23.*k32.*k34.*k43.*t2.*2.0+k14.*k21.*k23.*k34.*k41.*k43.*t2.*2.0)+k12.*k14.*k43+k12.*k23.*k41+k14.*k32.*k43-c.*k14.*k21.*k32+c.*k12.*k23.*k34+c.*k14.*k21.*k43+c.*k23.*k34.*k41+k21.*k23.*k34.*t2).*(-1.0./2.0))./(k12.*k14.*k23+c.*k14.*k21.*k23);
