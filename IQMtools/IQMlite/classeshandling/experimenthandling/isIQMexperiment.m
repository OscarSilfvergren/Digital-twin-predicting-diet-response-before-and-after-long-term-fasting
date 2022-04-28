function [output] = isIQMexperiment(input)
% isIQMexperiment: check if input argument is an IQMexperiment.

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>

output = strcmp(class(input),'IQMexperiment');
