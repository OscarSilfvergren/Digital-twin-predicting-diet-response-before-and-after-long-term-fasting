function [projectinfo] = parseMONOLIXprojectHeaderIQM(projectPath)
% Parses the project header information from the MONOLIX project and returns it.
% 
% [SYNTAX]
% [projectinfo] = parseMONOLIXprojectHeaderIQM(projectPath)
%
% [INPUT]
% projectPath:      Project to return the project header
%
% [OUTPUT]
% projectinfo:      Structure with information

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

% Check if project.mlxtran in project folder
if exist([projectPath '/project.mlxtran']),
    project = fileread([projectPath '/project.mlxtran']);
else
    error('project.mlxtran file could not be found.');
end

% Get the header
ixstart = strfind(project,'; ==PROJECT HEADER START===================================================');
ixend = strfind(project,  '; ==PROJECT HEADER END=====================================================');
if isempty(ixstart) || isempty(ixend),
    error('Project header could not be found in project.nmctl file.');
end
headertext = strtrim(project(ixstart+75:ixend-1));
headerterms = explodePCIQM(headertext,char(10));

% Construct output
projectinfo = [];
for k=1:length(headerterms),
    eval(['projectinfo.' strrep(strtrim(headerterms{k}(2:end)),'=','=explodePCIQM(') ','','',''['','']'');']);
end

% Add project.mlxtran to header
projectinfo.projectFile = 'project.mlxtran';
projectinfo.TOOL        = 'MONOLIX';

