function [boundaryPoints, metric_array_filtered, random_grid] = findBoundaryPoints(...
                          metric_array_curr,use_indices,min_points_per_bin,max_grid_res,ds_factor,reg_flag)

ignore_rand = false;                        
if ~exist('ds_factor')                        
    ds_factor = 100;
    ignore_rand = true;
end  
if ~exist('reg_factor')   
    reg_flag = false;
end    

% define mesh for edge sampling
metric_array_filtered = metric_array_curr(use_indices,:);
nPoints = size(metric_array_filtered,1);
nEdges = min(ceil(nPoints/min_points_per_bin),max_grid_res)+1;
if ~reg_flag
    grid1 = quantile(unique(metric_array_filtered(:,1)),linspace(0,1,nEdges));
    grid2 = quantile(unique(metric_array_filtered(:,2)),linspace(0,1,nEdges));
else
    grid1 = linspace(min(metric_array_filtered(:,1)),max(metric_array_filtered(:,1)),nEdges);
    grid2 = linspace(min(metric_array_filtered(:,2)),max(metric_array_filtered(:,2)),nEdges);
end    
% calculate number of random points to sample per grid patch

% assign obs to bins
[~,~,~,bin1,bin2] = histcounts2(metric_array_filtered(:,1),metric_array_filtered(:,2),grid1,grid2);   

% perform random grid-based sampling
if ~ignore_rand
    [~,~,ic] = unique([bin1 bin2],'rows');
    freq = accumarray(ic, 1);
    freq = freq(freq>0);
    n_rand_total = ceil(nPoints/ds_factor);
    random_grid = randsample(1:size(metric_array_filtered,1),n_rand_total,true,1./freq(ic));
else
    random_grid = [];
end    
% find edge values      
EdgeIndices1 = NaN(1,2*(nEdges-1));
% last_i = 1;
for b = 1:2:2*length(grid1)
    binIndices = find(bin2==ceil(b/2));
%     % random sampling
%     if length(binIndices)>1
%         rand_points = randsample(binIndices,min([n_rand_per_bin length(binIndices)]),false);
%     else
%         rand_points = binIndices;
%     end
%     random_grid(last_i+1:last_i+length(rand_points)) = rand_points;
%     last_i = last_i+length(rand_points);
    % gridded sampling
    if ~isempty(binIndices)
        [~, min_i] = min(metric_array_filtered(binIndices,1));
        EdgeIndices1(b) = binIndices(min_i);  
        [~, max_i] = max(metric_array_filtered(binIndices,1));
        EdgeIndices1(b+1) = binIndices(max_i);  
    end
end

EdgeIndices2 = NaN(1,2*(nEdges-1));
for b = 1:2:2*length(grid1)
    binIndices = find(bin1==ceil(b/2));
%     % random sampling
%     if length(binIndices)>1
%         rand_points = randsample(binIndices,min([n_rand_per_bin length(binIndices)]),false);
%     else
%         rand_points = binIndices;
%     end
%     random_grid(last_i+1:last_i+length(rand_points)) = rand_points;
%     last_i = last_i+length(rand_points);
    % gridded sampling
    if ~isempty(binIndices)
        [~, min_i] = min(metric_array_filtered(binIndices,2));
        EdgeIndices2(b) = binIndices(min_i);  
        [~, max_i] = max(metric_array_filtered(binIndices,2));
        EdgeIndices2(b+1) = binIndices(max_i);  
    end
end

boundaryPoints = [EdgeIndices1' ; EdgeIndices2'];  

% remove NaNs from random points
random_grid = random_grid(~isnan(random_grid))';