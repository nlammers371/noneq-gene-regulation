function [rates_out,error_flags] = ratesHM(rate_array,c_range,rate_bounds)

for k = 1:4
    eval(['k' num2str(k) ' = rate_array(:,k);'])
end
for r = 1:4
    eval(['r' num2str(r) ' = rate_array(:,r+4);'])
end

% k1Min = (((k4+r1+c_range.*r2).*r3.*r4) - c_range.*r1.*r2.*(c_range.*k2+r3+r4))./((k4+r1).*(c_range.*k2+r4));
% 
r1Half = (c_range.*k1.*k2+c_range.^2.*k2.*r2+c_range.*r2.*r3+k1.*r4+c_range.*r2.*r4+(-1).*r3.*r4) ...
  .^(-1).*(c_range.*k1.*k2.*k3+(-1).*c_range.*k1.*k2.*k4+(-1).*k1.*k3.*k4+c_range.* ...
  k2.*k3.*k4+c_range.^2.*k2.*k3.*r2+k3.*k4.*r3+c_range.*k3.*r2.*r3+(-1).*k1.* ...
  k4.*r4+k4.*r3.*r4+c_range.*r2.*r3.*r4);

k1Half = (c_range.*k2.*k3+(-1).*c_range.*k2.*k4+(-1).*k3.*k4+(-1).*c_range.*k2.*r1+(-1).*k4.* ...
  r4+(-1).*r1.*r4).^(-1).*((-1).*c_range.*k2.*k3.*k4+(-1).*c_range.^2.*k2.*k3.* ...
  r2+c_range.^2.*k2.*r1.*r2+(-1).*k3.*k4.*r3+(-1).*c_range.*k3.*r2.*r3+c_range.*r1.* ...
  r2.*r3+c_range.*r1.*r2.*r4+(-1).*k4.*r3.*r4+(-1).*r1.*r3.*r4+(-1).*c_range.* ...
  r2.*r3.*r4);

r2Half = c_range.^(-1).*(c_range.*k2.*k3+(-1).*c_range.*k2.*r1+k3.*r3+(-1).*r1.*r3+(-1).*r1.* ...
  r4+r3.*r4).^(-1).*((-1).*c_range.*k1.*k2.*k3+c_range.*k1.*k2.*k4+k1.*k3.*k4+( ...
  -1).*c_range.*k2.*k3.*k4+c_range.*k1.*k2.*r1+(-1).*k3.*k4.*r3+k1.*k4.*r4+k1.* ...
  r1.*r4+(-1).*k4.*r3.*r4+(-1).*r1.*r3.*r4);

k2Half = c_range.^(-1).*(k1.*k3+(-1).*k1.*k4+k3.*k4+(-1).*k1.*r1+c_range.*k3.*r2+(-1).* ...
  c_range.*r1.*r2).^(-1).*(k1.*k3.*k4+(-1).*k3.*k4.*r3+(-1).*c_range.*k3.*r2.* ...
  r3+c_range.*r1.*r2.*r3+k1.*k4.*r4+k1.*r1.*r4+c_range.*r1.*r2.*r4+(-1).*k4.* ...
  r3.*r4+(-1).*r1.*r3.*r4+(-1).*c_range.*r2.*r3.*r4);

r3Half = (k3.*k4+c_range.*k3.*r2+(-1).*c_range.*r1.*r2+k4.*r4+r1.*r4+c_range.*r2.*r4).^(-1).* ...
  ((-1).*c_range.*k1.*k2.*k3+c_range.*k1.*k2.*k4+k1.*k3.*k4+(-1).*c_range.*k2.*k3.*k4+ ...
  c_range.*k1.*k2.*r1+(-1).*c_range.^2.*k2.*k3.*r2+c_range.^2.*k2.*r1.*r2+k1.*k4.*r4+ ...
  k1.*r1.*r4+c_range.*r1.*r2.*r4);

k3Half = (c_range.*k1.*k2+(-1).*k1.*k4+c_range.*k2.*k4+c_range.^2.*k2.*r2+k4.*r3+c_range.*r2.*r3) ...
  .^(-1).*(c_range.*k1.*k2.*k4+c_range.*k1.*k2.*r1+c_range.^2.*k2.*r1.*r2+c_range.*r1.*r2.* ...
  r3+k1.*k4.*r4+k1.*r1.*r4+c_range.*r1.*r2.*r4+(-1).*k4.*r3.*r4+(-1).*r1.* ...
  r3.*r4+(-1).*c_range.*r2.*r3.*r4);

% if one or more are non-positive, select randomly
ref_vec = [1:3 5:7];
ref_array = repmat(ref_vec,size(rate_array,1),1);
% generate shuggling array
[~,col_array] = sort(rand(size(rate_array,1),6),2);
row_array = repmat((1:size(rate_array,1))',1,6);
lin_ind_array = sub2ind(size(rate_array),row_array,col_array);
% shuffle
ref_array = ref_array(lin_ind_array);
sol_array = [k1Half, k2Half, k3Half, r1Half, r2Half,r3Half];
sol_array = sol_array(lin_ind_array);
% filter
valid_bin_array = sol_array>=10^rate_bounds(1) & sol_array<=10^rate_bounds(2);
sol_array(~valid_bin_array) = NaN;
% ref_array(~valid_bin_array) = NaN;

% sample
[m_vals,m_ids] = nanmax(valid_bin_array,[],2);
lin_index_sub = sub2ind(size(ref_array),1:size(ref_array,1),m_ids');
col_index = ref_array(lin_index_sub);

% generate output
rate_vec_out = sol_array(lin_index_sub);
lin_index_out = sub2ind(size(rate_array),1:size(rate_array,1),col_index)';
error_flags = m_vals==0;
rates_out = rate_array;
rates_out(lin_index_out) = rate_vec_out;
% if length(pos_indices)>1
%   ind = randsample(pos_indices,1);
%   rate = sol_vec(ind);
% else



% (k1.*k2.*(k4+r3)+(k2+r1).*r3.*r4).*(k2.*k3.*k4+k4.*r1.*(k3+r2)+ ...
%   r1.*r2.*r3+k1.*(k4.*(k3+r2)+r2.*r3+k2.*(k3+k4+r3))+r1.*(k3+r2).* ...
%  r4+(r1+r2).*r3.*r4+k2.*(k3+r3).*r4).^(-1);
%    
 