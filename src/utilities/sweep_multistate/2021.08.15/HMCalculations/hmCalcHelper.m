function [delta,pON] = hmCalcHelper(RNum,activeStates)

    [V,D] = eig(RNum);
    [~,mi] = max(real(diag(D)));
    ss_vec_full = V(:,mi)./sum(V(:,mi));
    pON = sum(ss_vec_full(activeStates));
%     pOFF = sum(ss_vec_full(~activeStates));
    delta = (pON - 0.5);