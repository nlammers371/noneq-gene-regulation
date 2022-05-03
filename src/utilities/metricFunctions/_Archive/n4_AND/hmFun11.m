function out1 = hmFun11(c,k12,k14,k21,k23,k32,k34,k41,k43)
%HMFUN11
%    OUT1 = HMFUN11(C,K12,K14,K21,K23,K32,K34,K41,K43)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    02-Feb-2021 09:58:06

t2 = k12.^2;
t3 = k14.^2;
t4 = k21.^2;
t5 = k23.^2;
t6 = k32.^2;
t7 = k34.^2;
t8 = k41.^2;
t9 = k43.^2;
out1 = ((sqrt(t3.*t4.*t5+t3.*t4.*t6+t2.*t5.*t7+t3.*t4.*t9+t2.*t7.*t8+t4.*t6.*t9+t5.*t7.*t8+t6.*t7.*t8-k12.*k23.*t7.*t8.*2.0+k12.*k32.*t7.*t8.*2.0+k14.*k32.*t4.*t9.*2.0-k23.*k32.*t3.*t4.*2.0+k12.*k41.*t5.*t7.*2.0-k14.*k43.*t4.*t6.*2.0-k23.*k32.*t7.*t8.*2.0-k23.*k41.*t2.*t7.*2.0+k23.*k43.*t3.*t4.*2.0-k32.*k43.*t3.*t4.*2.0-k12.*k14.*k21.*k34.*t5.*2.0-k12.*k21.*k34.*k41.*t5.*4.0-k12.*k23.*k32.*k41.*t7.*2.0+k14.*k21.*k34.*k41.*t5.*2.0+k14.*k21.*k34.*k41.*t6.*2.0+k14.*k23.*k32.*k43.*t4.*2.0+k14.*k21.*k34.*k43.*t6.*4.0+k21.*k34.*k41.*k43.*t6.*2.0+k12.*k14.*k21.*k23.*k32.*k34.*2.0-k12.*k14.*k21.*k23.*k34.*k41.*2.0-k12.*k14.*k21.*k23.*k34.*k43.*2.0+k12.*k14.*k21.*k32.*k34.*k41.*2.0+k12.*k14.*k21.*k32.*k34.*k43.*4.0+k12.*k21.*k23.*k32.*k34.*k41.*4.0-k12.*k14.*k21.*k34.*k41.*k43.*2.0+k12.*k21.*k23.*k32.*k34.*k43.*2.0-k14.*k21.*k23.*k32.*k34.*k41.*4.0-k14.*k21.*k23.*k32.*k34.*k43.*4.0-k12.*k21.*k23.*k34.*k41.*k43.*4.0+k14.*k21.*k23.*k34.*k41.*k43.*2.0+k12.*k21.*k32.*k34.*k41.*k43.*2.0-k14.*k21.*k32.*k34.*k41.*k43.*2.0-k21.*k23.*k32.*k34.*k41.*k43.*2.0)+k14.*k21.*k23-k14.*k21.*k32+k12.*k23.*k34+k14.*k21.*k43-k12.*k34.*k41+k21.*k32.*k43+k23.*k34.*k41-k32.*k34.*k41).*(-1.0./2.0))./(c.*k21.*k23.*k34-c.*k21.*k32.*k34);
