function [result] = minIQM(varargin)
% minIQM: This function is used instead of the MATLAB "min" function, 
% allowing more than two scalar input arguments, each of type "double".
% 
% USAGE:
% ======
% [result] = minIQM(arg1,arg2,...,argn)   
%
% arg1...argn: scalar input arguments of type double.
%
% Output Arguments:
% =================
% result: min(arg1,arg2,arg3, ....)

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO IT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result = min([varargin{:}]);
return

