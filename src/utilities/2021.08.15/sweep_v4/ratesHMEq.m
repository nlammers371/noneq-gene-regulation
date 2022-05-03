function [rates_out,error_flags] = ratesHMEq(rate_array,c_range,rate_bounds)

for k = 1:4
    eval(['k' num2str(k) ' = rate_array(:,k);'])
end
for r = 1:4
    eval(['r' num2str(r) ' = rate_array(:,r+4);'])
end

% k1Min = (((k4+r1+c_range.*r2).*r3.*r4) - c_range.*r1.*r2.*(c_range.*k2+r3+r4))./((k4+r1).*(c_range.*k2+r4));
% 
f1a = (1/2).*k1.^(-1).*r1.^(-1).*(c_range.*k2+r4).^(-1).*((-1).*k1.*(c_range.*k2.*(( ...
  -1).*k3+k4)+k4.*(k3+r4))+(-1).*r1.*(c_range.^2.*k2.*r2+(-1).*r3.*r4+c_range.* ...
  r2.*(r3+r4))+(4.*k1.*r1.*(k4+c_range.*r2).*(c_range.*k2+r4).*(c_range.*k2.*k3+r3.*( ...
  k3+r4))+(c_range.*k1.*k2.*((-1).*k3+k4)+c_range.^2.*k2.*r1.*r2+(-1).*r1.*r3.* ...
  r4+k1.*k4.*(k3+r4)+c_range.*r1.*r2.*(r3+r4)).^2).^(1/2));

f1b = (-1/2).*k1.^(-1).*r1.^(-1).*(c_range.*k2+r4).^(-1).*((-1).*c_range.*k1.*k2.* ...
  k3+c_range.*k1.*k2.*k4+k1.*k3.*k4+c_range.^2.*k2.*r1.*r2+c_range.*r1.*r2.*r3+k1.* ...
  k4.*r4+c_range.*r1.*r2.*r4+(-1).*r1.*r3.*r4+(4.*k1.*r1.*(k4+c_range.*r2).*(c_range.* ...
  k2+r4).*(c_range.*k2.*k3+r3.*(k3+r4))+(c_range.*k1.*k2.*((-1).*k3+k4)+c_range.^2.* ...
  k2.*r1.*r2+(-1).*r1.*r3.*r4+k1.*k4.*(k3+r4)+c_range.*r1.*r2.*(r3+r4)) ...
  .^2).^(1/2));

f2a = (1/2).*c_range.^(-2).*k2.^(-1).*(k3+(-1).*r1).^(-1).*r2.^(-1).*(c_range.*((-1) ...
  .*k2.*k3.*k4+k1.*k2.*((-1).*k3+k4+r1)+r2.*((-1).*k3.*r3+(-1).*r3.* ...
  r4+r1.*(r3+r4)))+(-1).*(c_range.^2.*(4.*k2.*(k3+(-1).*r1).*r2.*(k1+(-1) ...
  .*r3).*(k3.*k4+(k4+r1).*r4)+(k2.*k3.*k4+k1.*k2.*(k3+(-1).*k4+(-1) ...
  .*r1)+r2.*(k3.*r3+r3.*r4+(-1).*r1.*(r3+r4))).^2)).^(1/2));

f2b = (1/2).*c_range.^(-2).*k2.^(-1).*(k3+(-1).*r1).^(-1).*r2.^(-1).*(c_range.*((-1) ...
  .*k2.*k3.*k4+k1.*k2.*((-1).*k3+k4+r1)+r2.*((-1).*k3.*r3+(-1).*r3.* ...
  r4+r1.*(r3+r4)))+(c_range.^2.*(4.*k2.*(k3+(-1).*r1).*r2.*(k1+(-1).*r3).* ...
  (k3.*k4+(k4+r1).*r4)+(k2.*k3.*k4+k1.*k2.*(k3+(-1).*k4+(-1).*r1)+ ...
  r2.*(k3.*r3+r3.*r4+(-1).*r1.*(r3+r4))).^2)).^(1/2));

f3a = (-1/2).*k3.^(-1).*(k4+c_range.*r2).^(-1).*r3.^(-1).*((-1).*k1.*k3.*k4+ ...
  c_range.*k2.*k3.*(k1+k4)+c_range.^2.*k2.*k3.*r2+k4.*r3.*r4+r1.*r3.*r4+c_range.*r2.* ...
  r3.*((-1).*r1+r4)+(4.*k3.*(k4+c_range.*r2).*(k1.*(k4+r1)+c_range.*r1.*r2).* ...
  r3.*(c_range.*k2+r4)+((-1).*k1.*k3.*k4+c_range.*k2.*k3.*(k1+k4)+c_range.^2.*k2.*k3.* ...
  r2+(k4+r1).*r3.*r4+c_range.*r2.*r3.*((-1).*r1+r4)).^2).^(1/2));

f3b = (-1/2).*k3.^(-1).*(k4+c_range.*r2).^(-1).*r3.^(-1).*((-1).*k1.*k3.*k4+ ...
  c_range.*k2.*k3.*(k1+k4)+c_range.^2.*k2.*k3.*r2+k4.*r3.*r4+r1.*r3.*r4+c_range.*r2.* ...
  r3.*((-1).*r1+r4)+(-1).*(4.*k3.*(k4+c_range.*r2).*(k1.*(k4+r1)+c_range.*r1.* ...
  r2).*r3.*(c_range.*k2+r4)+((-1).*k1.*k3.*k4+c_range.*k2.*k3.*(k1+k4)+c_range.^2.* ...
  k2.*k3.*r2+(k4+r1).*r3.*r4+c_range.*r2.*r3.*((-1).*r1+r4)).^2).^(1/2));

% if one or more are non-positive, select randomly
ref_vec = [repelem(1:3,2)];
ref_array = repmat(ref_vec,size(rate_array,1),1);
ref_rates1 = repelem(rate_array(:,1:3),1,2);
ref_rates2 = repelem(rate_array(:,5:7),1,2);

% generate shuffling array
[~,col_array] = sort(rand(size(rate_array,1),6),2);
row_array = repmat((1:size(rate_array,1))',1,6);
lin_ind_array = sub2ind(size(rate_array),row_array,col_array);

% shuffle and check for imaginary components
ref_array = ref_array(lin_ind_array);
ref_rates = cat(3,ref_rates1(lin_ind_array),ref_rates2(lin_ind_array));

real_flag_array = [imag(f1a)==0, imag(f1b)==0, imag(f2a)==0, imag(f2b)==0, imag(f3a)==0, imag(f3b)==0];
sol_array = [real(f1a), real(f1b), real(f2a), real(f2b), real(f3a), real(f3b)];
sol_array = sol_array(lin_ind_array);
sol_array = repmat(sol_array,1,1,2) .* ref_rates;

% filter
valid_bin_array = nanmin(real_flag_array & sol_array>=10^rate_bounds(1) & sol_array<=10^rate_bounds(2),[],3);
sol1 = sol_array(:,:,1);
sol2 = sol_array(:,:,2);

sol1(~valid_bin_array) = NaN;
sol2(~valid_bin_array) = NaN;

% sample
[m_vals,m_ids] = nanmax(nanmin(valid_bin_array,[],3),[],2);
lin_index_sub = sub2ind(size(ref_array),1:size(ref_array,1),m_ids');
col_index = ref_array(lin_index_sub);

% generate output
rate_vec_out = [sol1(lin_index_sub)' ; sol2(lin_index_sub)'];%repmat(sol_array(lin_index_sub)',2,1);
lin_index_out = [sub2ind(size(rate_array),1:size(rate_array,1),col_index)' ;
                 sub2ind(size(rate_array),1:size(rate_array,1),col_index+4)'] ;
error_flags = m_vals==0;
rates_out = rate_array;
rates_out(lin_index_out) = rate_vec_out;