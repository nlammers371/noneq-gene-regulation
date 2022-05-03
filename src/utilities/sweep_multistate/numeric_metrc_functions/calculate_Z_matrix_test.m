function Z_num = calculate_Z_matrix_test(Q_num,ss_num,num_prec)

    % make sure the inputs are correct
%     assert(round(sum(ss_num),4)==1,'SS vec is not normalized');
%     assert(all(round(sum(Q_num),4)==0),'Q matrix is of improper form');
    if round(sum(ss_num),num_prec)==1 && all(round(sum(Q_num),num_prec)==0)
        % uses eq 55 from Whitt 1992 to calculate Z (the "fundamental matrix"
        % of Q
        ss_num = reshape(ss_num,1,[]);

        % make PI helper matrix
        PI = repmat(ss_num,size(Q_num,1),1);

        % calculate Z
        Z_num = (PI - Q_num')\eye(size(Q_num,1)) - PI;
    else
        Z_num = NaN(size(Q_num));
    end