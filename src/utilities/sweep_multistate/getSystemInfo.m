function simInfo = getSystemInfo
    
    % note that this assumes that 1 and only 1 networkInfo mat file is in
    % the working path
    load('networkInfo.mat','networkInfo')
    
    % initialize info struct
    simInfo = struct;
    
    % transfer fields to simInfo
    fnames = fieldnames(networkInfo);
    for f = 1:length(fnames)
        simInfo.(fnames{f}) = networkInfo.(fnames{f});
    end
    