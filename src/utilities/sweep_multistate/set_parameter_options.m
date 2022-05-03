function simInfo = set_parameter_options(simInfo,simType)

% extract basic info
crs = simInfo.crs;
% determine if we are considering off-target interference
% ns_ind = strfind(simInfo.functionPath,'_ns');
% n_ns = str2double(simInfo.functionPath(ns_ind+3:ns_ind+4));
specOnlyFlag = ~any(simInfo.n_total_bound~=simInfo.n_right_bound);

if ~isfield(simInfo,'sweepVarStrings')
    simInfo.sweepVarStrings = cellstr(string(simInfo.sweepVarList));
end    
% always ensure this equality
simInfo.defaultValues(1) = crs;

constraintPairCell = {};

% get param indices

if simInfo.nStates > 4 && ~specOnlyFlag
    if isfield(simInfo,'a_index')
        index_vec = [simInfo.cr_index simInfo.cw_index simInfo.a_index];
    else
        index_vec = [simInfo.cr_index simInfo.cw_index simInfo.b_index];
    end
    cw = simInfo.cw;
    specFactor = simInfo.specFactor;
    simInfo.defaultValues(index_vec(1:3)) = [crs cw specFactor];
else  
    index_vec = [simInfo.cr_index];
    simInfo.defaultValues(index_vec(1)) = crs;
end
% 
% if strcmp(simType,'all_specific') && simInfo.nStates > 4
%     simInfo.defaultValues(index_vec(1:3)) = [crs realmin 1]; % note that we set specificity factor to "1" and cw=cr
%     simInfo.cr_index = 1;  
%     
% elseif strcmp(simType,'PIC') && simInfo.nStates > 4 % standard setup for >4 state system (which for now by definition includes non-specific factors
%     % set binding rates of PIC to be constant
%     kap_index = find(strcmp(simInfo.sweepVarStrings,'kap'));
%     kam_index = find(strcmp(simInfo.sweepVarStrings,'kam'));
%     constraintPairCell(1) = {[kap_index kam_index]}; 
%     
%     % set all cooperativity terms to 1
%     wap_indices = contains(simInfo.sweepVarStrings,'wap');
%     simInfo.defaultValues(wap_indices) = 1;
%     simInfo.sweepFlags(wap_indices) = false;
%     
%     % set binding rates of TFs to be equivalent, regardless of PIC status
%     kpi_index = find(strcmp(simInfo.sweepVarStrings,'kpi'));
%     kpa_index = find(strcmp(simInfo.sweepVarStrings,'kpa'));
%     constraintPairCell(2) = {[kpi_index kpa_index]};     
%     
% elseif strcmp(simType,'Nucleosome') && simInfo.nStates > 4 % standard setup for >4 state system (which for now by definition includes non-specific factors
%   
%     % set binding rates of nucleosome to be constant
%     kip_index = find(strcmp(simInfo.sweepVarStrings,'kip'));
%     kim_index = find(strcmp(simInfo.sweepVarStrings,'kim'));
%     constraintPairCell(1) = {[kip_index kim_index]}; 
%     
%     % set all cooperativity terms to 1
%     wip_indices = contains(simInfo.sweepVarStrings,'wip');
%     simInfo.defaultValues(wip_indices) = 1;
%     simInfo.sweepFlags(wip_indices) = false;
%     
%     % set unbinding rates of TFs to be equivalent, regardless of PIC status
%     kmi_index = find(strcmp(simInfo.sweepVarStrings,'kmi'));
%     kma_index = find(strcmp(simInfo.sweepVarStrings,'kma'));
%     constraintPairCell(2) = {[kmi_index kma_index]}; 
%     
%     % check for hgher order stuff
%     for n = 2:10
%         wmi_index = find(strcmp(simInfo.sweepVarStrings,['wmi' num2str(n)]));
%         wma_index = find(strcmp(simInfo.sweepVarStrings,['wma' num2str(n)]));
%         if ~isempty(wma_index)
%             constraintPairCell(end+1) = {[wmi_index wma_index]}; 
%         end
%     end
% end
simInfo.simType = simType;
simInfo.constraintPairCell = constraintPairCell;
simInfo.specOnlyFlag = specOnlyFlag;
simInfo.numCalcFlag = contains(simInfo.functionPath,'numeric');

% deal with "collapse" options
% if strcmp(simInfo.modelCollapseOption,'generic_act')
%     simInfo.map_to_indices = simInfo.to_indices_generic_act;
%     simInfo.map_from_indices = simInfo.from_indices_generic_act;
% elseif strcmp(simInfo.modelCollapseOption,'generic_gen')
%     simInfo.map_to_indices = simInfo.to_indices_generic_gen;
%     simInfo.map_from_indices = simInfo.from_indices_generic_gen;    
% else    
%     simInfo.map_to_indices = [];
%     simInfo.map_from_indices = [];    
% end    