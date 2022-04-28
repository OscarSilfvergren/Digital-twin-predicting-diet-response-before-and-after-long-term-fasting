function [output] = isIQMmodel(model)
% isIQMmodel: check if input argument is an IQMmodel.

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>

output = strcmp(class(model),'IQMmodel');
