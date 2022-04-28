function [ output ] = isMONOLIXprojectIQM( projectPath )
% Checks if given project path is a Monolix project. This is checked by 
% requiring a project.mlxtran file in this folder.
% 
% [SYNTAX]
% [output] = isMONOLIXprojectIQM( projectPath )
%
% [INPUT]
% projectPath:  PAth to check if Monolix project
%
% [OUTPUT]
% output = 0: No MONOLIX project
% output = 1: MONOLIX project

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

testfile = fullfile(projectPath,'project.mlxtran');
if ~exist(testfile),
    output = 0;
else
    output = 1;
end


