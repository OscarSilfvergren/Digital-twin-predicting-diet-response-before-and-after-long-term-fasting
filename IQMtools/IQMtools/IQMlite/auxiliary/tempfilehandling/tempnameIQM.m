function [tmpname] = tempnameIQM( input_args )
% tempnameIQM does the same as tempname in MATLAB but uses the IQM Tools own
% TEMP folder.

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>

[dummy, filename] = fileparts(tempname);
tmpname = [tempdirIQM,filename];
