function var_num = calculate_var_num(Q_num,ss_num,activeStates,Z_num,num_prec)


    % make sure the inputs are correct
%     assert(round(sum(ss_num),4)==1,'SS vec is not normalized');
%     assert(all(round(sum(Q_num),4)==0),'Q matrix is of improper form');
    if round(sum(ss_num),num_prec)==1 && all(round(sum(Q_num),num_prec)==0)
        % now caculate the variance (see eq 12 from above from Whitt 1992)
        ss_num = reshape(ss_num,1,[]);
        % generate emissions helper vec
        f_vec = zeros(1,size(Q_num,1)); % initiation rate for each state
        f_vec(activeStates) = 1; % assume all active states produce at the same rate

        % use matrix form
        ZF = Z_num*f_vec';
        var_num = 2*(ss_num.*f_vec)*ZF;
        
%         %%
%         r_vec = zeros(1,size(Q_num,1));
%         r_vec(activeStates) = 1/3;
%         ZFR = Z_num*r_vec';
%         var_num_r = 2*(ss_num.*r_vec)*ZFR;
%         %%
    else
        var_num = NaN;
    end

%     varSS = 0;
%     for i = 1:size(Q_num,1)
%         for j = 1:size(Q_num,2)
%             varSS = varSS + f_vec(i)*ss_num(i)*Z_num(i,j)*f_vec(j);
%         end
%     end
%     varSS = 2*varSS
