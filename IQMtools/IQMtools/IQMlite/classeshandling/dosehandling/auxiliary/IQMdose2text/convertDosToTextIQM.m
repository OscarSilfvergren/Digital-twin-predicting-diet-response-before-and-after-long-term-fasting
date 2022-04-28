function [dosTextStructure] = convertDosToTextIQM(dos)
% convertDosToTextIQM: Converts an IQMdosing object to a structure
% containing the different parts of the text description of the dosing
% scheme.

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>

% Initialize variables
dosTextStructure = [];
% Get IQMstructure
ds = IQMstruct(dos);
% Parse structure into the dosTextStructure description

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dosTextStructure.name = ds.name;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dosTextStructure.notes = ds.notes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dosTextStructure.inputs = {};
allInputText = {};
for k=1:length(ds.inputs),
    inputText = '';
    input = ds.inputs(k);
    % Write the limiter (name)
    inputText = sprintf('********** %s\n',input.name);
    % Type
    inputText = sprintf('%stype:         %s\n',inputText,input.type);
    % Time
    if length(input.time) == 1,
        inputText = sprintf('%stime:         %g\n',inputText,input.time);
    else
        timeText = sprintf('%g, ',input.time);
        inputText = sprintf('%stime:         [%s]\n',inputText,timeText(1:end-2));
    end
    % Tlag (write only if non-empty or non-zero)
    if ~isempty(input.Tlag),
        if ischar(input.Tlag),
            if isempty(input.TlagNotes),
                inputText = sprintf('%sTlag:         %s\n',inputText,input.Tlag);
            else
                inputText = sprintf('%sTlag:         %s    %% %s\n',inputText,input.Tlag,input.TlagNotes);
            end
        elseif isnumeric(input.Tlag),
            if input.Tlag~=0,
                if isempty(input.TlagNotes),
                    inputText = sprintf('%sTlag:         %g\n',inputText,input.Tlag);
                else
                    inputText = sprintf('%sTlag:         %g    %% %s\n',inputText,input.Tlag,input.TlagNotes);
                end
            end
        end
    end
    % D
    if length(input.D) == 1,
        inputText = sprintf('%sD:            %g\n',inputText,input.D);
    else
        doseText = sprintf('%g, ',input.D);
        inputText = sprintf('%sD:            [%s]\n',inputText,doseText(1:end-2));
    end
    % parameters
    if ~isempty(input.parameters),
        for k2=1:length(input.parameters),
            parametersValueText = sprintf('%g, ',input.parameters(k2).value);
            if isempty(input.parameters(k2).notes),
                inputText = sprintf('%s%s:%s[%s]\n',inputText,input.parameters(k2).name,char(32*ones(1,14-length(input.parameters(k2).name)-1)),parametersValueText(1:end-2));
            else
                inputText = sprintf('%s%s:%s[%s]    %% %s\n',inputText,input.parameters(k2).name,char(32*ones(1,14-length(input.parameters(k2).name)-1)),parametersValueText(1:end-2),input.parameters(k2).notes);
            end
        end
    end
    % notes
    if ~isempty(input.notes),
        inputText = sprintf('%snotes:        %s\n',inputText,input.notes);
    end
    % save text
    dosTextStructure.inputs{k} = inputText;
end

