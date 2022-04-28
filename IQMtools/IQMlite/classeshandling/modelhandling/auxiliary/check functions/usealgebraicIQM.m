function [output] = usealgebraicIQM(model)
% usealgebraicIQM: checks if an IQMmodel contains algebraic rules.
%
% Output Arguments:
% =================
% output: =1 if algebraic rules are used, =0 if not

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>

if ~isIQMmodel(model),
    error('Given model is not an IQMmodel.');
end

% default setting
output = 0; % no algebraic rule present

% get model structure
ms = struct(model);

% check algebraic rules
if ~isempty(ms.algebraic),
    output = 1;
end
