function [] = IQMcreateTEXTmodel(name,varargin)
% IQMcreateTEXTmodel: Creates a new TEXT model, saves the corresponding file
% in the current directory and opens it in the editor, if opening is
% desired.
%
% USAGE:
% ======
% IQMcreateTEXTmodel(name)
% IQMcreateTEXTmodel(name,openFlag)
%
% name: filename of the model and model name
% openFlag: decides if the created model is automatically opened in the
%           editor or not. (=0: do not open, =1: do open)
%
% DEFAULT VALUES:
% ===============
% openFlag: 1 (do open)

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLE VARIABE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
openFlag = 1;
if nargin == 2,
    openFlag = varargin{1};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE NEW MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = IQMmodel();
ms = struct(model);
ms.name = name;
model = IQMmodel(ms);
IQMcreateTEXTfile(model,name);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPEN IF DESIRED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if openFlag,
    edit(strcat(name,'.txt'));
end

return
