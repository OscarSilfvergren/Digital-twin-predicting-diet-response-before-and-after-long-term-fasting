function [output] = isvariableIQM(model,name)
% isvariableIQM: checks if "name" is a variable in the provided model.
% This function works only for IQMmodels. The check is of
% course case sensitive
%
% Output Arguments:
% =================
% output: =1 if "name" is a variable, =0 if not

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>


if ~strcmp(class(model),'IQMmodel'),
    error('Given model is not an IQMmodel.');
end

% get all variables of the model
allVariables = IQMvariables(model);

% check if "name" is a variable.
output = 0;
for k = 1:length(allVariables),
    if strcmp(strtrim(name),allVariables{k}),
        output = 1;
        break;
    end
end
return