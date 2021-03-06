function a12Sym = a12FunFromMathematica_v2

syms cr cw b kmi kma kpi kpa kap kam kip kim positive

a12Sym = (1/2).*(b.*cr.^2.*kap.^2.*kma.*kpa.*kpi+cr.*cw.*kap.^2.*kma.*kpa.* ...
  kpi+b.*cr.*cw.*kap.^2.*kma.*kpa.*kpi+cw.^2.*kap.^2.*kma.*kpa.*kpi+( ...
  -1).*b.*cr.^2.*kap.*kip.*kma.*kpa.*kpi+(-1).*cr.*cw.*kap.*kip.*kma.* ...
  kpa.*kpi+(-1).*b.*cr.*cw.*kap.*kip.*kma.*kpa.*kpi+(-1).*cw.^2.* ...
  kap.*kip.*kma.*kpa.*kpi+b.*cr.^2.*kap.*kip.*kmi.*kpa.*kpi+cr.*cw.* ...
  kap.*kip.*kmi.*kpa.*kpi+b.*cr.*cw.*kap.*kip.*kmi.*kpa.*kpi+cw.^2.* ...
  kap.*kip.*kmi.*kpa.*kpi+(-1).*b.*cr.^2.*kip.^2.*kmi.*kpa.*kpi+(-1) ...
  .*cr.*cw.*kip.^2.*kmi.*kpa.*kpi+(-1).*b.*cr.*cw.*kip.^2.*kmi.*kpa.* ...
  kpi+(-1).*cw.^2.*kip.^2.*kmi.*kpa.*kpi+b.^2.*cr.^2.*kap.*kma.*kmi.* ...
  kpa.*kpi+2.*b.*cr.*cw.*kap.*kma.*kmi.*kpa.*kpi+cw.^2.*kap.*kma.* ...
  kmi.*kpa.*kpi+(-1).*b.^2.*cr.^2.*kip.*kma.*kmi.*kpa.*kpi+(-2).*b.* ...
  cr.*cw.*kip.*kma.*kmi.*kpa.*kpi+(-1).*cw.^2.*kip.*kma.*kmi.*kpa.* ...
  kpi).^(-1).*((-1).*b.*cr.*kam.*kap.^2.*kma.*kpa+(-1).*cw.*kam.* ...
  kap.^2.*kma.*kpa+b.*cr.*kam.*kap.*kip.*kma.*kpa+cw.*kam.*kap.*kip.* ...
  kma.*kpa+(-1).*b.*cr.*kam.*kap.*kip.*kmi.*kpa+(-1).*cw.*kam.*kap.* ...
  kip.*kmi.*kpa+b.*cr.*kam.*kip.^2.*kmi.*kpa+cw.*kam.*kip.^2.*kmi.* ...
  kpa+(-1).*b.*cr.*kam.*kap.*kma.*kmi.*kpa+(-1).*b.^2.*cr.*kam.*kap.* ...
  kma.*kmi.*kpa+(-1).*cw.*kam.*kap.*kma.*kmi.*kpa+(-1).*b.*cw.*kam.* ...
  kap.*kma.*kmi.*kpa+b.^2.*cr.*kam.*kip.*kma.*kmi.*kpa+cw.*kam.*kip.* ...
  kma.*kmi.*kpa+b.*cr.*kap.*kip.*kma.*kmi.*kpa+b.*cw.*kap.*kip.*kma.* ...
  kmi.*kpa+(-1).*b.*cr.*kam.*kip.*kmi.^2.*kpa+(-1).*b.*cw.*kam.*kip.* ...
  kmi.^2.*kpa+b.*cr.*kip.^2.*kmi.^2.*kpa+b.*cw.*kip.^2.*kmi.^2.*kpa+( ...
  -1).*b.^2.*cr.*kam.*kma.*kmi.^2.*kpa+(-1).*b.*cw.*kam.*kma.* ...
  kmi.^2.*kpa+b.^2.*cr.*kip.*kma.*kmi.^2.*kpa+b.*cw.*kip.*kma.* ...
  kmi.^2.*kpa+(-1).*b.*cr.*kap.^2.*kim.*kma.*kpi+(-1).*cw.*kap.^2.* ...
  kim.*kma.*kpi+b.*cr.*kap.*kim.*kip.*kma.*kpi+cw.*kap.*kim.*kip.* ...
  kma.*kpi+(-1).*b.*cr.*kap.^2.*kma.^2.*kpi+(-1).*b.*cw.*kap.^2.* ...
  kma.^2.*kpi+b.*cr.*kap.*kim.*kma.^2.*kpi+b.*cw.*kap.*kim.*kma.^2.* ...
  kpi+(-1).*b.*cr.*kap.*kim.*kip.*kmi.*kpi+(-1).*cw.*kap.*kim.*kip.* ...
  kmi.*kpi+b.*cr.*kim.*kip.^2.*kmi.*kpi+cw.*kim.*kip.^2.*kmi.*kpi+( ...
  -1).*b.^2.*cr.*kap.*kim.*kma.*kmi.*kpi+(-1).*cw.*kap.*kim.*kma.* ...
  kmi.*kpi+(-1).*b.*cr.*kap.*kip.*kma.*kmi.*kpi+(-1).*b.*cw.*kap.* ...
  kip.*kma.*kmi.*kpi+b.*cr.*kim.*kip.*kma.*kmi.*kpi+b.^2.*cr.*kim.* ...
  kip.*kma.*kmi.*kpi+cw.*kim.*kip.*kma.*kmi.*kpi+b.*cw.*kim.*kip.* ...
  kma.*kmi.*kpi+(-1).*b.^2.*cr.*kap.*kma.^2.*kmi.*kpi+(-1).*b.*cw.* ...
  kap.*kma.^2.*kmi.*kpi+b.^2.*cr.*kim.*kma.^2.*kmi.*kpi+b.*cw.*kim.* ...
  kma.^2.*kmi.*kpi+((b.*cr.*kam.*kap.^2.*kma.*kpa+cw.*kam.*kap.^2.* ...
  kma.*kpa+(-1).*b.*cr.*kam.*kap.*kip.*kma.*kpa+(-1).*cw.*kam.*kap.* ...
  kip.*kma.*kpa+b.*cr.*kam.*kap.*kip.*kmi.*kpa+cw.*kam.*kap.*kip.* ...
  kmi.*kpa+(-1).*b.*cr.*kam.*kip.^2.*kmi.*kpa+(-1).*cw.*kam.*kip.^2.* ...
  kmi.*kpa+b.*cr.*kam.*kap.*kma.*kmi.*kpa+b.^2.*cr.*kam.*kap.*kma.* ...
  kmi.*kpa+cw.*kam.*kap.*kma.*kmi.*kpa+b.*cw.*kam.*kap.*kma.*kmi.* ...
  kpa+(-1).*b.^2.*cr.*kam.*kip.*kma.*kmi.*kpa+(-1).*cw.*kam.*kip.* ...
  kma.*kmi.*kpa+(-1).*b.*cr.*kap.*kip.*kma.*kmi.*kpa+(-1).*b.*cw.* ...
  kap.*kip.*kma.*kmi.*kpa+b.*cr.*kam.*kip.*kmi.^2.*kpa+b.*cw.*kam.* ...
  kip.*kmi.^2.*kpa+(-1).*b.*cr.*kip.^2.*kmi.^2.*kpa+(-1).*b.*cw.* ...
  kip.^2.*kmi.^2.*kpa+b.^2.*cr.*kam.*kma.*kmi.^2.*kpa+b.*cw.*kam.* ...
  kma.*kmi.^2.*kpa+(-1).*b.^2.*cr.*kip.*kma.*kmi.^2.*kpa+(-1).*b.* ...
  cw.*kip.*kma.*kmi.^2.*kpa+b.*cr.*kap.^2.*kim.*kma.*kpi+cw.*kap.^2.* ...
  kim.*kma.*kpi+(-1).*b.*cr.*kap.*kim.*kip.*kma.*kpi+(-1).*cw.*kap.* ...
  kim.*kip.*kma.*kpi+b.*cr.*kap.^2.*kma.^2.*kpi+b.*cw.*kap.^2.* ...
  kma.^2.*kpi+(-1).*b.*cr.*kap.*kim.*kma.^2.*kpi+(-1).*b.*cw.*kap.* ...
  kim.*kma.^2.*kpi+b.*cr.*kap.*kim.*kip.*kmi.*kpi+cw.*kap.*kim.*kip.* ...
  kmi.*kpi+(-1).*b.*cr.*kim.*kip.^2.*kmi.*kpi+(-1).*cw.*kim.*kip.^2.* ...
  kmi.*kpi+b.^2.*cr.*kap.*kim.*kma.*kmi.*kpi+cw.*kap.*kim.*kma.*kmi.* ...
  kpi+b.*cr.*kap.*kip.*kma.*kmi.*kpi+b.*cw.*kap.*kip.*kma.*kmi.*kpi+( ...
  -1).*b.*cr.*kim.*kip.*kma.*kmi.*kpi+(-1).*b.^2.*cr.*kim.*kip.*kma.* ...
  kmi.*kpi+(-1).*cw.*kim.*kip.*kma.*kmi.*kpi+(-1).*b.*cw.*kim.*kip.* ...
  kma.*kmi.*kpi+b.^2.*cr.*kap.*kma.^2.*kmi.*kpi+b.*cw.*kap.*kma.^2.* ...
  kmi.*kpi+(-1).*b.^2.*cr.*kim.*kma.^2.*kmi.*kpi+(-1).*b.*cw.*kim.* ...
  kma.^2.*kmi.*kpi).^2+(-4).*(b.*kam.*kap.^2.*kma.^2+(-1).*b.* ...
  kap.^2.*kim.*kma.^2+2.*b.*kam.*kap.*kip.*kma.*kmi+(-2).*b.*kap.* ...
  kim.*kip.*kma.*kmi+b.*kam.*kap.*kma.^2.*kmi+b.^2.*kam.*kap.* ...
  kma.^2.*kmi+(-1).*b.*kap.*kim.*kma.^2.*kmi+(-1).*b.^2.*kap.*kim.* ...
  kma.^2.*kmi+b.*kam.*kip.^2.*kmi.^2+(-1).*b.*kim.*kip.^2.*kmi.^2+ ...
  b.*kam.*kip.*kma.*kmi.^2+b.^2.*kam.*kip.*kma.*kmi.^2+(-1).*b.* ...
  kim.*kip.*kma.*kmi.^2+(-1).*b.^2.*kim.*kip.*kma.*kmi.^2+b.^2.* ...
  kam.*kma.^2.*kmi.^2+(-1).*b.^2.*kim.*kma.^2.*kmi.^2).*(b.*cr.^2.* ...
  kap.^2.*kma.*kpa.*kpi+cr.*cw.*kap.^2.*kma.*kpa.*kpi+b.*cr.*cw.* ...
  kap.^2.*kma.*kpa.*kpi+cw.^2.*kap.^2.*kma.*kpa.*kpi+(-1).*b.*cr.^2.* ...
  kap.*kip.*kma.*kpa.*kpi+(-1).*cr.*cw.*kap.*kip.*kma.*kpa.*kpi+(-1) ...
  .*b.*cr.*cw.*kap.*kip.*kma.*kpa.*kpi+(-1).*cw.^2.*kap.*kip.*kma.* ...
  kpa.*kpi+b.*cr.^2.*kap.*kip.*kmi.*kpa.*kpi+cr.*cw.*kap.*kip.*kmi.* ...
  kpa.*kpi+b.*cr.*cw.*kap.*kip.*kmi.*kpa.*kpi+cw.^2.*kap.*kip.*kmi.* ...
  kpa.*kpi+(-1).*b.*cr.^2.*kip.^2.*kmi.*kpa.*kpi+(-1).*cr.*cw.* ...
  kip.^2.*kmi.*kpa.*kpi+(-1).*b.*cr.*cw.*kip.^2.*kmi.*kpa.*kpi+(-1).* ...
  cw.^2.*kip.^2.*kmi.*kpa.*kpi+b.^2.*cr.^2.*kap.*kma.*kmi.*kpa.*kpi+ ...
  2.*b.*cr.*cw.*kap.*kma.*kmi.*kpa.*kpi+cw.^2.*kap.*kma.*kmi.*kpa.* ...
  kpi+(-1).*b.^2.*cr.^2.*kip.*kma.*kmi.*kpa.*kpi+(-2).*b.*cr.*cw.* ...
  kip.*kma.*kmi.*kpa.*kpi+(-1).*cw.^2.*kip.*kma.*kmi.*kpa.*kpi)).^( ...
  1/2));
