function probSuccess =  probSuccess(rLow,rHigh,sigmaLow,sigmaHigh,T)

r_critical = (sigmaHigh.*rLow + sigmaLow.*rHigh) ./ (sigmaLow+sigmaHigh);

probSuccess = ...
    1 - .5*(.5*(1 - erf((r_critical-rLow)./sqrt(2)./ sigmaLow .* sqrt(T))) + ...
    .5*(1-erf((rHigh-r_critical) ./sqrt(2)./ sigmaHigh .* sqrt(T))));