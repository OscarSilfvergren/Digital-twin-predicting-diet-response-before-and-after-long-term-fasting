function [IQMstructure] = IQMstruct(dos)
% IQMstruct: This function returns the internal data structure
% of an IQMdosing object
%
% USAGE:
% ======
% [IQMstructure] = IQMstruct(dos) 
%
% dos: IQMdosing object 
%
% Output Arguments:
% =================
% IQMstructure: internal data structure of the IQMdosing object

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VERY SIMPLE FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IQMstructure = struct(dos);