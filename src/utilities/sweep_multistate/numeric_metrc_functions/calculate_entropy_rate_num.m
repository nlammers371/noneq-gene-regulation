function phi_num = calculate_entropy_rate_num(Q_num,ss_num,num_prec)


    % make sure the inputs are correct
%     assert(round(sum(ss_num),4)==1,'SS vec is not normalized');
%     assert(all(round(sum(Q_num),4)==0),'Q_num matrix is of improper form');
    phi_num = NaN;
    if round(sum(ss_num),num_prec)==1 && all(round(sum(Q_num),num_prec)==0)
        % now caculate the variance (see eq 12 from Whitt 1992)
        ss_num = reshape(ss_num,1,[]);

        % use matrix form for caculations
        QR = Q_num ./ Q_num';
        QR(isnan(QR)) = 1;    
        QM = Q_num .* log(QR);

        phi_num = sum(ss_num*QM');
    end
%     phi1 = 0;
%     for i = 1:size(Q_num,1)
%         for j = 1:size(Q_num,1)
%             if i~=j && Q_num(i,j)~=0
%                 phi1 = phi1 + ss_num(i)*Q_num(j,i)*log(Q_num(j,i)/Q_num(i,j));
%             end
%         end
%     end
        

    