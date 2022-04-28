function [project] = IQMupdatemeasurement(project,experimentindex,measurement,varargin)
% IQMupdatemeasurement: update or add a measurement in an experiment of a project.  
%
% USAGE:
% ======
% [project] = IQMupdatemeasurement(project,experimentindex,measurement)
% [project] = IQMupdatemeasurement(project,experimentindex,measurement,measurementindex)
%
% project:          IQMprojectSB object
% experimentindex:  index of the experiment to add the measurement to
% measurement:      IQMmeasurement which to update or add
% measurementindex: index of the measurement to be updated. If omitted the measurment is
%                   added to the experiment as last measurement.
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
if ~isnumeric(experimentindex),
    error('Second input argument is not an experiment index.');
end
if ~isIQMmeasurement(measurement),
    error('Third input argument is not an IQMmeasurement.');
end
projectstruct = IQMstruct(project);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 3,
    measurementindex = length(projectstruct.experiments(experimentindex).measurements)+1;
elseif nargin == 4,
    measurementindex = varargin{1};
    if measurementindex < 1 || measurementindex > length(projectstruct.experiments(experimentindex).measurements),
        error('''measurementindex'' out of bounds.');
    end
else
    error('Incorrect number of input arguments.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adding/Updating the project
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectstruct.experiments(experimentindex).measurements{measurementindex} = measurement;
project = IQMprojectSB(projectstruct);
return
