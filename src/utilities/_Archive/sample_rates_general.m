function new_rate_vec = sample_rates_general(current_rate_vec,rate_bounds,sigma)          
    % initialize vector to store new rates    
    new_rate_vec = NaN(size(current_rate_vec));
    for i = 1:numel(current_rate_vec)
%         pd = makedist('LogNormal',log(current_rate_vec(i)),sigma);        
%         pd = truncate(pd,rate_bounds(1),rate_bounds(2));
        pd = makedist('LogNormal',0,sigma);        
        pd = truncate(pd,rate_bounds(1)/current_rate_vec(i),rate_bounds(2)/current_rate_vec(i));
        new_rate_vec(i) = current_rate_vec(i)*random(pd,1);
    end    
end  