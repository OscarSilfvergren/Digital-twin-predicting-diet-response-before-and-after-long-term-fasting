function [output] = isIQMdosing(input)
% isIQMdosing: check if input argument is an IQMdosing object.

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>

output = strcmp(class(input),'IQMdosing');
