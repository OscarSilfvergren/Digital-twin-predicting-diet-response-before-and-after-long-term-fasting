%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD PARAMETER INDICES AND NOMINAL VALUES TO WORKESTIMATION STRUCTURE
% Could be different indices for different experiments ... so its done here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [workestimation] = addparameterindices(workestimation,parameters,parameterslocal)
for k=1:length(workestimation),
    [allparameters,nominalparamvalues] = IQMparameters(workestimation(k).IQMmodel);
    paramindices = getnamevectorindices(allparameters,parameters.names);
    paramindiceslocal = getnamevectorindices(allparameters,parameterslocal.names);
    workestimation(k).paramindices = paramindices;
    workestimation(k).paramindiceslocal = paramindiceslocal;
    workestimation(k).paramnominal = nominalparamvalues;
end
return

