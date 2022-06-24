function RSymWrong = RSymFunWrong(cw,b,ki,ka,km,kp,wip,wap,wma,wpa)
%RSYMFUNWRONG
%    RSYMWRONG = RSYMFUNWRONG(CW,B,KI,KA,KM,KP,WIP,WAP,WMA,WPA)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    23-Jun-2022 11:56:49

t2 = b.*km;
t3 = cw.*kp;
t4 = ka.*wap;
t5 = ki.*wip;
t6 = t2.*wma;
t7 = t3.*wpa;
RSymWrong = reshape([-ka-t3,t3,0.0,ka,t2,-t2-t4,t4,0.0,0.0,t5,-t5-t6,t6,ki,0.0,t7,-ki-t7],[4,4]);
