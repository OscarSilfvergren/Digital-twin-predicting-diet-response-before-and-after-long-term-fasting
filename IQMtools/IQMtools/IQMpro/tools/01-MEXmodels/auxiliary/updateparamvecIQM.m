function [newparamvector] = updateparamvecIQM(paramnames, oldparamvector, pnames, pvalues)
% updateparamvecIQM: The mex simulation functions take a full parameter value
% vector as an input. In order to construct such a vector from an existing 
% parameter vector while changing some but possibly not all parameter
% values this function can be used. 
% 
% USAGE:
% ======
% [newparamvector] = updateparamvecIQM(paramnames, oldparamvector, pnames, pvalues)
% 
% paramnames: Cell-array of all model parameter names (in the order as they
%   appear in the model)
% oldparamvector: Full parameter vector before the change
% pnames: Cell-array of parameter names for which to change the
%   parametervalues 
% pvalues: Vector with parameter values for the parameters to be changed
%   (order defined by pnames)
% 
% OUTPUT:
% =======
% newparamvector: Constructed new parameter vector with changed parameters

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

% Convert to cell-array if necessary
if ischar(pnames),
    pnames = {pnames};
end
% Parameter input arguments need to have same length ...
if length(pnames) ~= length(pvalues),
    error('Different numbers of parameters and parametervalues to be changed.');
end

% Update the paramvector with given values
newparamvector = oldparamvector;
for k = 1:length(pnames),
    index = strmatchIQM(pnames{k}, paramnames, 'exact');
    if isempty(index),
        disp(sprintf('Parameter name ''%s'' is no parameter in the model.',pnames{k}));
    end
    newparamvector(index) = pvalues(k);
end
return