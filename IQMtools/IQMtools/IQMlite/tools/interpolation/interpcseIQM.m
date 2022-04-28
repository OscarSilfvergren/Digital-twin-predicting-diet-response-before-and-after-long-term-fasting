% interpcseIQM: Cubic spline interpolation with endpoints
% This is a MEX function, so this .m file only contains the documentation.
%
% USAGE:
% ======
% yy = interpcseIQM(x,y,xx)
% yy = interpcseIQM(x,y,xx,e1,e2)
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
