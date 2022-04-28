function [output] = isIQMprojectSB(project)
% isIQMprojectSB: check if input argument is an IQMprojectSB

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

output = strcmp(class(project),'IQMprojectSB');
