function [IQMstructure] = IQMstruct(iqm)
% IQMstruct: This function returns the internal data structure
% of an IQMexperiment object
%
% USAGE:
% ======
% [IQMstructure] = IQMstruct(experiment) 
%
% experiment: IQMexperiment object 
%
% Output Arguments:
% =================
% IQMstructure: internal data structure of the IQMexperiment object

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VERY SIMPLE FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IQMstructure = struct(iqm);