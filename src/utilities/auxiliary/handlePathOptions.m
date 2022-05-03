function writePath = handlePathOptions(subFolderName)

    % figure out where we are
    currentDir = pwd;
    % slashes = sort([strfind(currentDir,'\') strfind(currentDir,'/')]);
    src_ind = strfind(currentDir,'src');
    utilities_path = [currentDir(1:src_ind+2) filesep 'utilities' filesep];
    
    writePath = [utilities_path 'metricFunctions' filesep subFolderName filesep];
    mkdir(writePath);

    metric_fun_path = [utilities_path 'metricFunctions' filesep];
    metric_script_path = [currentDir(1:src_ind+2) filesep 'metric_calculation_scripts' filesep];

    addpath(genpath(utilities_path))
    addpath(genpath(metric_script_path))
    rmpath(genpath(metric_fun_path));
    addpath(writePath);