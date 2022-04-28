function [yy] = interpcsexIQM(x,y,xx,varargin)
% interpcsexIQM: Cubic spline interpolation with endpoints
% This is an interface function only to interpcseIQM. Only in the C-code
% simulation of models there is a difference. interpcseIQM does only
% evaluate the spline coefficients once in the first call. Subsequent calls
% only evaluate the spline function. This leads to a considerable speed-up.
% On the other hand, interpcsexIQM in the C-code simulation does evaluate
% the spline coefficients each time it is called. => time dependent splines
% can be implemented. However, for MATLAB simulation there is no
% difference.
%
% USAGE:
% ======
% yy = interpcsexIQM(x,y,xx)
% yy = interpcsexIQM(x,y,xx,e1,e2)
% 
% x:     x-values 
% y:     y-values
% xx:    x-values at which to evaluate the spline function (allow multiple)
% e1,e2: endpoint derivatives (if specified, both need to be given)
% yy:    interpolated value
%
% DEFAULT VALUES:
% ===============
% e1, e2:   =0 (no endpoint slopes defined)
%
% Output Arguments:
% =================
% yy:    scalar or vector of interpolated values.

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>


if nargin < 3 || nargin > 5,
    error('interpcsexIQM: incorrect number of input arguments.');
elseif nargin == 3,
    yy = interpcseIQM(x,y,xx);
elseif nargin == 5,
    yy = interpcseIQM(x,y,xx,varargin{1},varargin{2});
end    
