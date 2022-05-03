function [new_rate_vec, Kd] = sample_rate_params(unbinding_edges,...
                    binding_edges, current_rate_vec, rate_bounds,c,varargin)
    iter = 0;
    sigma = .25;
    exit_flag = false;
    while ~exit_flag
        iter = iter + 1;
        % initialize vector to store new rates    
        prop_rate_vec = NaN(size(current_rate_vec));
        % naively draw new values subject to constraints
        for i = 1:numel(prop_rate_vec)
            pd = makedist('LogNormal',log(current_rate_vec(i)),sigma);
            ub = rate_bounds(2);
            if ismember(i,binding_edges)
                ub = ub / c;
            end
            pd = truncate(pd,rate_bounds(1),ub);
            prop_rate_vec(i) = random(pd,1);
        end
        % enforce consistent Kd
        prop_rate_vec(unbinding_edges(2)) = prop_rate_vec(binding_edges(2)) * ...
            prop_rate_vec(unbinding_edges(1)) / prop_rate_vec(binding_edges(1));
        % enforce detailed balance
        options = find(~ismember(1:8,[unbinding_edges binding_edges]));
        norm_edge = randsample(options,1);
        if norm_edge > 4
            norm_value = prod(prop_rate_vec(1:4))/prod(prop_rate_vec(5:8)) * prop_rate_vec(norm_edge);
        else
            norm_value = prod(prop_rate_vec(5:8))/prod(prop_rate_vec(1:4)) * prop_rate_vec(norm_edge);
        end
        prop_rate_vec(norm_edge) = norm_value;        
        % check for violations
        boundary_violations = any(prop_rate_vec>rate_bounds(2) | prop_rate_vec<rate_bounds(1));   
        if ~boundary_violations
            new_rate_vec = prop_rate_vec;
            Kd = new_rate_vec(binding_edges(1))/new_rate_vec(unbinding_edges(1));
            exit_flag = true;
        elseif iter > 10
            new_rate_vec = NaN(size(prop_rate_vec));
            Kd = new_rate_vec(binding_edges(1))/new_rate_vec(unbinding_edges(1));
            exit_flag = true;
            warning('unable to find suitable sample')
        end
    end  