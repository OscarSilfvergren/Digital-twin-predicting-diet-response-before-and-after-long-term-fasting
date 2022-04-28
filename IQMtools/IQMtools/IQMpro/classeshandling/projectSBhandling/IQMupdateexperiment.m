function [project] = IQMupdateexperiment(project,experiment,varargin)
% IQMupdateexperiment: update or add a experiment in a project.  
%
% USAGE:
% ======
% [project] = IQMupdateexperiment(project,experiment)
% [project] = IQMupdateexperiment(project,experiment,experimentindex)        
%
% project:         IQMprojectSB object
% experiment:      IQMexperiment which to update or add
% experimentindex: index of the experiment to be updated. If omitted the experiment is
%                  added to the project as last experiment.
%
% Output Arguments:
% =================
% project: updated project

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check input argument
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isIQMprojectSB(project),
    error('First input argument is not an IQMprojectSB.');
end
if ~isIQMexperiment(experiment),
    error('Second input argument is not an IQMexperiment.');
end
projectstruct = IQMstruct(project);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 2,
    experimentindex = length(projectstruct.experiments)+1;
elseif nargin == 3,
    experimentindex = varargin{1};
    if experimentindex < 1 || experimentindex > length(projectstruct.experiments),
        error('''experimentindex'' out of bounds.');
    end
else
    error('Incorrect number of input arguments.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adding/Updating the project
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectstruct.experiments(experimentindex).experiment = experiment;
if isempty(projectstruct.experiments(experimentindex).name),
    % use experiments name as name for the experiment in the project
    % (otherwise empty if added experiments)
    x = struct(experiment);
    projectstruct.experiments(experimentindex).name = x.name;
end
project = IQMprojectSB(projectstruct);
return
