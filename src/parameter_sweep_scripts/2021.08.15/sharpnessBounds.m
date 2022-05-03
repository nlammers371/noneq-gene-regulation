

beta = 100;
cr = 1;
cw = 1e4;
mu = cr/cw;
f_eq = beta;
f_vec = logspace(log10(beta),log10(beta^2));
% sharpnessBound = f_vec./(4*beta).*((beta-1).*(f_vec+beta))./(f_vec.^2-1);

% sharpnessBound2 = (1/4).*f_vec.*(1+f_vec).^(-1).*(1+(1+(-1).*(beta.*cr+cw).*(beta.*cr+(-1).*cw.*f_vec).* ...
%                 (beta.^2.*cr+(-1).*cw.*f_vec).^(-1)).^(-1));
%               
sharpnessSimp = ((-1 + beta).* (beta + f_vec))./(beta .* (-1 + f_vec));

sharpnessBound = (1/4).*mu.*((-1)+beta).*beta.^(-1).*((-1)+f_vec).^(-1).*f_vec.*(beta+f_vec).*(1+mu.* ...
                f_vec).^(-1);
              
sharpnessBoundApp = (1/4).*mu.*(f_vec).^(-1).*f_vec.*(beta+f_vec).*(1+mu.* ...
                f_vec).^(-1);              
              
close all
figure;
hold on
plot(sharpnessSimp,log10(f_vec./f_eq))
% plot(log10(f_vec./f_eq),sharpnessBoundApp)
% plot(log10(f_vec./f_eq),sharpnessBound)