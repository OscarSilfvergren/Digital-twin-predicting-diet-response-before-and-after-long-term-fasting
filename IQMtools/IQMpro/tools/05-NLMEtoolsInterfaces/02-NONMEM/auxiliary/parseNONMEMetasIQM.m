function [ dataeta, OMEGA, OMEGAnames ] = parseNONMEMetasIQM( projectPath )
% Parses a NONMEM project and returns the ETAs.
% 
% [SYNTAX]
% [ dataeta, OMEGA, OMEGAnames ] = parseNONMEMetasIQM( projectPath )
%
% [INPUT]
% projectPath:      Project to return the ETA information
%
% [OUTPUT]
% dataeta:          ETA data
% OMEGA:            OMEGA data
% OMEGAnames:       Names of omegas

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

% Construct RESULTS path
resultsPath = [projectPath '/RESULTS'];
    
% Check the projectPath
if ~exist(resultsPath),
    error(sprintf('The provided project path "%s" does not point to a valid NONMEM project.\nPlease make sure a "RESULTS" folder is in the provided path.',projectPath));
end

% Check that indiv_eta.txt is present in the RESULTS folder
indiv_eta_file = [resultsPath '/project.eta'];
if ~exist(indiv_eta_file)
    error('The "project.eta" file does not exist in the RESULTS folder.');
end

% Determine random effect estimates for shrinkage determination
x = parseNONMEMresultsIQM(projectPath);
y = sampleNONMEMpopulationParametersIQM(x,0,1);
OMEGA       = y.randomEffects.values;
OMEGAnames  = y.randomEffects.names;

% Load eta file
indiv_eta   = IQMloadNONCSVdataset([resultsPath '/project.eta'],1);

% Get eta modes
dataeta = table();
for k=1:length(OMEGAnames),
    dataeta.(OMEGAnames{k}) = indiv_eta.(['ETA_' OMEGAnames{k}]);
end

