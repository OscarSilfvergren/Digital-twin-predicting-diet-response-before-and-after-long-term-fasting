function [structure] = IQMstruct(project)
% IQMstruct: This function returns the internal data structure
% of an IQMprojectSB
%
% USAGE:
% ======
% [structure] = IQMstruct(IQMprojectSB) 
%
% IQMprojectSB: IQMprojectSB 
%
% Output Arguments:
% =================
% structure: internal data structure of the IQMprojectSB

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VERY SIMPLE FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
structure = struct(project);