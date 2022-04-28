function [output] = isIQMmeasurement(input)
% isIQMmeasurement: check if input argument is an IQMmeasurement.

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>

output = strcmp(class(input),'IQMmeasurement');
