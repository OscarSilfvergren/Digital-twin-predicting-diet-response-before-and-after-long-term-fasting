function [output] = isparameterIQM(model,name)
% isparameterIQM: checks if "name" is a parameter in the provided model.
% This function works both for IQMmodels and ODE files. The check is of
% course case sensitive
%
% Output Arguments:
% =================
% output: =1 if "name" is a parameter, =0 if not

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>


% get all parameters of the model
allParameters = IQMparameters(model);

% check if "name" is a parameter.
output = 0;
for k = 1:length(allParameters),
    if strcmp(strtrim(name),allParameters{k}),
        output = 1;
        break;
    end
end
return