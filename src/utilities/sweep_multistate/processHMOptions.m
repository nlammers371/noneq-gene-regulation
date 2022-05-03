function simInfo = processHMOptions(simInfo)
    
    if simInfo.half_max_flag
        if ismember(simInfo.nStates,[4 6]) && ~simInfo.numCalcFlag
            loadPath = [simInfo.functionPath 'HMFunctions' filesep];
            load([loadPath 'hmStruct.mat'],'hmStruct');
            simInfo.hmStruct = hmStruct;
        end
    end
   
    