function [baseNum, activeStatesBase, n_wrong_bound_init,n_right_bound_init,...
        activity_vec_init,transitionInfoInit,permittedConnectionsInit]...
            = generateBasicSixStateNetwork

    % define some basic parameters
    activeStatesBase = [3 4 5];
    baseNum = 6;

    % zero out forbidden transitions
    [fromRef,toRef] = meshgrid(1:baseNum,1:baseNum);
    diffRef = abs(fromRef-toRef);
    toRefHalved = toRef<=baseNum/2;
    fromRefHalved = fromRef<=baseNum/2;
    permittedConnectionsRaw= (diffRef==1 & toRefHalved==fromRefHalved) | diffRef==baseNum/2;

    % permute these connections to follow a more intuitive labeling scheme
    indexCurr = 1:baseNum;
    indexAdjusted = circshift(indexCurr,floor(baseNum/4));
    indexAdjusted = [indexAdjusted(1:baseNum/2) fliplr(indexAdjusted(baseNum/2+1:end))];
    [~,si] = sort(indexAdjusted);
    permittedConnectionsInit = permittedConnectionsRaw(si,si);

    % generate an array with binding info
    transitionInfoInit = zeros(size(permittedConnectionsInit));

    % specific binding/unbinding
    spec_pairs = {[1,2],[4,3]};
    for s = 1:length(spec_pairs)
        ind1 = spec_pairs{s}(1);
        ind2 = spec_pairs{s}(2);
        % update
        transitionInfoInit(ind2,ind1) = 1;
        transitionInfoInit(ind1,ind2) = -1;
    end

    % non-specific binding/unbinding
    non_spec_pairs = {[1,6],[4,5]};
    for s = 1:length(non_spec_pairs)
        ind1 = non_spec_pairs{s}(1);
        ind2 = non_spec_pairs{s}(2);
        % update
        transitionInfoInit(ind2,ind1) = 2;
        transitionInfoInit(ind1,ind2) = -2;
    end

    % locus activity fluctuations
    locus_pairs = {[1,4],[6,5],[2,3]};
    for s = 1:length(locus_pairs)
        ind1 = locus_pairs{s}(1);
        ind2 = locus_pairs{s}(2);
        % update
        transitionInfoInit(ind2,ind1) = 3;
        transitionInfoInit(ind1,ind2) = -3;
    end

    % generate array that indicates activity state
    activity_vec_init = false(1,baseNum);
    activity_vec_init(activeStatesBase) = 1;

    % generate flags indicating numbers of right and wrong factors bound
    n_right_bound_init = zeros(size(activity_vec_init));
    n_wrong_bound_init = n_right_bound_init;
    n_right_bound_init([2 3]) = 1;
    n_wrong_bound_init([5 6]) = 1;
    
end    