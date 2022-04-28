function [result] = iscomplex (x)
% Checks if number x is a complex number.
%
% USAGE:
% ======
% [result] = iscomplex (x)

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>

result = ~isreal(x);