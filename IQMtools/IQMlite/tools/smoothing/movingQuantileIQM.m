function [y]=movingQuantileIQM(x,range,q)
% movingQuantileIQM will compute moving quantiles over a range, defined by the user.
%
% Usage: y = movingQuantileIQM(x,range,q)
% where 
%      x:       is the input vector (or matrix) to be smoothed. 
%      range: 	percentage of number of points in x to be used to average over
%      q:       quantileIQM to be calculated
%      y:       is output vector of same length as x
%
% Note:if x is a matrix then the smoothing will be done 'vertically'.
%
% This function is a wrapper for the movingAverageIQM function.
% 
% Examples:
%   x=randn(300,1);
%   plot(x,'g.'); 
%   hold on;
%   plot(movingQuantileIQM(x,30,0.2),'r--');

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>

[y] = movingAverageIQM(x,range,@(x)quantileIQM(x,q));



