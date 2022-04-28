function [datastruct] = createdatastructIQMplotIQM(varargin)
% createdatastructIQMplotIQM: Creates a datastruct for the IQMplot function
%
% USAGE:
% ======
% [datastruct] = createdatastructIQMplotIQM(time,data)
% [datastruct] = createdatastructIQMplotIQM(time,data,names)
% [datastruct] = createdatastructIQMplotIQM(time,data,names,name)
% [datastruct] = createdatastructIQMplotIQM(time,data,names,legendtext,name)
% [datastruct] = createdatastructIQMplotIQM(time,data,names,legendtext,marker,name)
% [datastruct] = createdatastructIQMplotIQM(time,data,names,errorindices,minvalues,maxvalues,legendtext,marker,name)
%
% time: column vector with time information
% data: matrix with data where each row corresponds to one time point and
%       each column to a different variable
% names: cell-array with the names of the data variables
%
% name: name for the datastruct
%
% legendtext: cell-array of same length as names with text to be used for
%             the legend.
% marker: marker and line style for plot
% errorindices: indices of the data for which errorbounds are available
% minvalues: error bounds for data ... to be shown by error bars
% maxvalues: error bounds for data ... to be shown by error bars
%  
% DEFAULT data:
% ===============
% names: the plotted variables obtain the name 'x1', 'x2', ...
% legendtext: same as names
% marker: '-'
% min/maxvalues: no errorbars shown
% name: 'unnamed'
%
% Output Arguments:
% =================
% datastruct: structure that can be displayed by IQMplot   (>> IQMplot(datastruct))

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initializing the datastruct structure and setting default values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datastruct = [];
datastruct.time = varargin{1};
datastruct.data = varargin{2};
datastruct.name = 'unnamed';
datastruct.dataNames = {};
for k = 1:size(datastruct.data,2);
    datastruct.dataNames{end+1} = sprintf('x%d',k);
end
datastruct.errorindices = [];
datastruct.minvalues = [];
datastruct.maxvalues = [];
datastruct.marker = '-';

if nargin == 2,
    % do nothing ... all done already
    datastruct.legendtext = datastruct.dataNames;
elseif nargin == 3,
    datastruct.dataNames = varargin{3};
    datastruct.legendtext = datastruct.dataNames;
elseif nargin == 4,
    datastruct.dataNames = varargin{3};
    datastruct.name = varargin{4};
    datastruct.legendtext = datastruct.dataNames;
elseif nargin == 5,
    datastruct.dataNames = varargin{3};
    datastruct.legendtext = varargin{4};
    datastruct.name = varargin{5};
elseif nargin == 6,
    datastruct.dataNames = varargin{3};
    datastruct.legendtext = varargin{4};
    datastruct.marker = varargin{5};
    datastruct.name = varargin{6};
elseif nargin == 9,
    datastruct.dataNames = varargin{3};
    datastruct.errorindices = varargin{4};
    datastruct.minvalues = varargin{5};
    datastruct.maxvalues = varargin{6};
    datastruct.legendtext = varargin{7};
    datastruct.marker = varargin{8};
    datastruct.name = varargin{9};
else
    error('Wrong number of input arguments.');
end
if isempty(datastruct.legendtext),
    datastruct.legendtext = datastruct.dataNames;
end


% Check data consistency
if size(datastruct.time,1) ~= size(datastruct.data,1),
    error('Different number of time points and time points in data.');
end
if length(datastruct.dataNames) ~= size(datastruct.data,2),
    error('Different number of variable data and variable names.');
end
