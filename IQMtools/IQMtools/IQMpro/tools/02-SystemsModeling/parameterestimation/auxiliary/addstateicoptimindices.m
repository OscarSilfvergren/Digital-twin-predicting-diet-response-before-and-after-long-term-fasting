%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD INDICES OF THE STATES FOR WHICH TO OPTIMIZE THE INITIAL CONDITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [workestimation] = addstateicoptimindices(workestimation,initialconditions)
for k=1:length(workestimation),
    stateindicesicoptim = getnamevectorindices(IQMstates(workestimation(k).model),initialconditions.names);
    workestimation(k).stateindicesicoptim = stateindicesicoptim;
end
