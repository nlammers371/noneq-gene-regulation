function RSymFull = RSymFunFull(cr,cw,a,ki,ka,kp,km,wap,wip,wma,wpa,wmp,wa,wi)
%RSYMFUNFULL
%    RSYMFULL = RSYMFUNFULL(CR,CW,A,KI,KA,KP,KM,WAP,WIP,WMA,WPA,WMP,WA,WI)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    17-Feb-2022 12:26:36

RSymFull = reshape([-ka-cr.*kp-cw.*kp,cr.*kp,0.0,ka,0.0,cw.*kp,km,-km-ka.*wap,ka.*wap,0.0,0.0,0.0,0.0,ki.*wip,-ki.*wip-km.*wma,km.*wma,0.0,0.0,ki,0.0,cr.*kp.*wpa,-ki-cr.*kp.*wpa-cw.*kp.*wpa,cw.*kp.*wpa,0.0,0.0,0.0,0.0,a.*km.*wma,-ki.*wip-a.*km.*wma,ki.*wip,a.*km,0.0,0.0,0.0,ka.*wap,-a.*km-ka.*wap],[6,6]);
