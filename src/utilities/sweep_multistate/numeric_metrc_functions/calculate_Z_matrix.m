function Z_num = calculate_Z_matrix(Q_num,ss_num,num_prec)

    % make sure the inputs are correct

    if round(sum(ss_num),num_prec)==1 && all(round(sum(Q_num),num_prec)==0)
        % uses eq 55 from Whitt 1992 to calculate Z (the "fundamental
        % matrix")
        % of Q
        ss_num = reshape(ss_num,1,[]);

        % make PI helper matrix
        PI = repmat(ss_num,size(Q_num,1),1);

        % calculate Z
        Z_num = inv(PI - Q_num') - PI;
    else
        Z_num = NaN(size(Q_num));
    end