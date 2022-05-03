function [Tau_ON_num,Tau_OFF_num,cycle_time] = calculate_tau_num(Q_num,ss_num,onStateFilter,num_prec)

    warning('error', 'MATLAB:nearlySingularMatrix');
    % make sure the inputs are correct
%     assert(round(sum(ss_num),4)==1,'SS vec is not normalized');
%     assert(all(round(sum(Q_num),4)==0),'Q matrix is of improper form');
    if round(sum(ss_num),num_prec)==1 && all(round(sum(Q_num),num_prec)==0)
        % now caculate the variance (see eq 12 from above from Whitt 1992)
        ss_num = reshape(ss_num,1,[]);

        %%%%%%%%%% ON RATE (off state dwell time) %%%%%%%%%%%%%%%%%%%%%
        offStateFilter = ~onStateFilter;

        % create truncated transition rate matrix
        QNumTruncON = zeros(sum(offStateFilter)+1,sum(offStateFilter)+1);
        QNumTruncON(1:end-1,1:end-1) = Q_num(offStateFilter,offStateFilter);
        QNumTruncON(end,1:end-1) = sum(Q_num(~offStateFilter,offStateFilter));
        QNumTruncON(1:end-1,end) = sum(Q_num(offStateFilter,~offStateFilter),2);
        QNumTruncON(end,end) = -sum(QNumTruncON(:,end));

        % calculate SS
    %     ss_vec_on = zeros(1,size(QNumTruncON,1));
    %     ss_vec_on(1:end-1) = ss_num(offStateFilter);
    %     ss_vec_on(end) = sum(ss_num(onStateFilter));
        ss_vec_on = calculate_ss_num(QNumTruncON,num_prec);

        % calculate Z
        Z_mat_on = calculate_Z_matrix(QNumTruncON,ss_vec_on,num_prec);

        % calculate ON times
        pON = ss_vec_on(end);
        ETON_vec = (Z_mat_on(end,end) - Z_mat_on(1:end-1,end))/pON;

        % calculate flux-weighted mean
        inFluxVecOFF = Q_num(offStateFilter,~offStateFilter) * ss_num(~offStateFilter)';
        Tau_OFF_num = (inFluxVecOFF'*ETON_vec) / sum(inFluxVecOFF);       


        %%%%%%%%%% OFF RATE (on state dwell time) %%%%%%%%%%%%%%%%%%%%%

        % create truncated transition rate matrix
        QNumTruncOFF = zeros(sum(onStateFilter)+1,sum(onStateFilter)+1);
        QNumTruncOFF(1:end-1,1:end-1) = Q_num(onStateFilter,onStateFilter);
        QNumTruncOFF(end,1:end-1) = sum(Q_num(offStateFilter,onStateFilter));
        QNumTruncOFF(1:end-1,end) = sum(Q_num(onStateFilter,offStateFilter),2);
        QNumTruncOFF(end,end) = -sum(QNumTruncOFF(:,end));

        % calculate SS
    %     ss_vec_off = zeros(1,size(QNumTruncOFF,1));
    %     ss_vec_off(1:end-1) = ss_num(onStateFilter);
    %     ss_vec_off(end) = sum(ss_num(offStateFilter));
        ss_vec_off = calculate_ss_num(QNumTruncOFF,3);

        % calculate Z
        Z_mat_off = calculate_Z_matrix(QNumTruncOFF,ss_vec_off,num_prec);

        % calculate ON times
        pOFF = ss_vec_off(end);
        ETOFF_vec = (Z_mat_off(end,end) - Z_mat_off(1:end-1,end))/pOFF;

        % calculate flux-weighted mean
        inFluxVecON = Q_num(onStateFilter,offStateFilter) * ss_num(offStateFilter)';
        Tau_ON_num = (inFluxVecON'*ETOFF_vec) / sum(inFluxVecON);

       % calculate cycle time
       cycle_time = Tau_OFF_num+Tau_ON_num;
    else
        Tau_ON_num = NaN;
        Tau_OFF_num = NaN;
        cycle_time = NaN;
    end