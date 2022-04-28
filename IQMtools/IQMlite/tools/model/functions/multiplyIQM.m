function [result] = multiplyIQM(varargin)
% multiplyIQM: Implementing the multiply MathML function
% 
% USAGE:
% ======
% [result] = multiplyIQM(arg1,arg2,...,argn)   

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>

result = 1;
for k = 1:nargin,
    result = result * varargin{k};
end
return

