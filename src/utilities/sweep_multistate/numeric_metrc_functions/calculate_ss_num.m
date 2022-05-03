function ss_num = calculate_ss_num(Q_num,num_prec)

    [V,D] = eig(Q_num);
    [mv,mi] = max(real(diag(D)));
    % check form
%     assert(round(mv,4)==0);
    ss_num = NaN(size(Q_num,1),1);
    if round(mv,num_prec)==0
        % normalize
        ss_num = V(:,mi)./sum(V(:,mi));
%         if nanmin(ss_num)<10^-num_prec
%             ss_num = NaN(size(ss_num));
%         end
%         ss_num = exp(log(V(:,mi))-logsumexp(log(V(:,mi)),1));
    end