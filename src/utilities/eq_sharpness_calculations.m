function [p2or4,Sharpness] = eq_sharpness_calculations(ec,em)


Sharpness = exp(1).^(ec+em).*((-1)+exp(1).^(ec+em)).*(exp(1).^ec+exp(1).^em+ ...
  2.*exp(1).^(ec+em)).^(-2);

p2or4 = 2.*(2+exp(1).^((-1).*ec)+exp(1).^((-1).*em)).^(-1);