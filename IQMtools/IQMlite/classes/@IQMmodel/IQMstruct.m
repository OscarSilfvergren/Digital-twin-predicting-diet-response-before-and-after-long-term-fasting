function [IQMstructure] = IQMstruct(iqm)
% IQMstruct: This function returns the internal data structure
% of an IQMmodel
%
% USAGE:
% ======
% [IQMstructure] = IQMstruct(IQMmodel) 
%
% IQMmodel: IQMmodel 
%
% Output Arguments:
% =================
% IQMstructure: internal data structure of the IQMmodel

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VERY SIMPLE FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IQMstructure = struct(iqm);