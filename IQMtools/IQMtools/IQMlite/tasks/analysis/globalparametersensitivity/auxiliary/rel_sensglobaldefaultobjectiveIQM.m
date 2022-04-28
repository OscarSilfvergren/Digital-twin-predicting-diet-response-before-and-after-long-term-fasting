function [rawsensitivities] = rel_sensglobaldefaultobjectiveIQM(MEXmodel,timevector,paramNames,paramValues,ICs,nomsimdata,stateindices,variableindices,reactionindices,integratoroptions)
% rel_sensglobaldefaultobjectiveIQM: Default global sensitivity analysis
% objective function. This function determines the sqrt of the sum of
% RELATIVE squared errors (SSSE), between simulations of the nominal model
% and the perturbed model.
%
% The output argument contains both the SSSEs of the Single components (SSSSE)
% (states, variables, and reactions) and a Combined version (CSSSE), which 
% appears as last element in the output vector.

global Ncount
Ncount = Ncount + 1;

% Simulate the model for given perturbations in the parameters
if isIQMproPresent(),
    pertsimdata = IQMPsimulate(MEXmodel,timevector,ICs,paramNames,paramValues,integratoroptions);
else
    modelSim = IQMparameters(MEXmodel,paramNames,paramValues);
    pertsimdata = IQMsimulate(modelSim,integratoroptions.method,timevector,ICs,integratoroptions);
end
% Get the nominal and perturbed elements for comparison
nomstatevalues = nomsimdata.statevalues(:,stateindices);
pertstatevalues = pertsimdata.statevalues(:,stateindices);
nomvariablevalues = nomsimdata.variablevalues(:,variableindices);
pertvariablevalues = pertsimdata.variablevalues(:,variableindices);
nomreactionvalues = nomsimdata.reactionvalues(:,reactionindices);
pertreactionvalues = pertsimdata.reactionvalues(:,reactionindices);
allnomvalues = [nomstatevalues nomvariablevalues nomreactionvalues];
allpertvalues = [pertstatevalues pertvariablevalues pertreactionvalues];
% get the relative differences
% For this we need to handle zero nominal values. We do that by replacing
% these zero nominal values by "1". The rational is that zero nominal
% values happen 
%   a) due to initial conditions (but then also the pert values are 0)
%   b) an output is always zero (but then so are the pert values)
%   c) well and otherwise it will become an approximation ...
allnomvaluesNonZero = allnomvalues;
allnomvaluesNonZero(find(allnomvalues==0)) = 1;
relativevalues = (allpertvalues-allnomvalues)./allnomvaluesNonZero;
% get the squared and summed differences
SSSSE = sqrt(sum((relativevalues).^2));
CSSSE = sqrt(sum(sum((relativevalues).^2))); 
% return the data
rawsensitivities = [SSSSE CSSSE]; 

