function [modelTextBCStructure] = convertModelToTextBCIQM(iqm)
% convertModelToTextBCIQM: Converts an IQMmodel to a structure containing the 
% different parts of the biochemical oriented text description of the
% model.

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>


% first check if the complete stoichiometric matrix can be determined.
% otherwise a biochemical representation is not fully possible and some
% states still need to be defined by differential equations.
% determine the names of that states that need to be described by ODEs
rawFlag = 1;
silentFlag = 1;
[N,stateNamesBC] = IQMstoichiometry(iqm,rawFlag,silentFlag);
stateNames = IQMstates(iqm);
stateNamesODE = setdiff(stateNames,stateNamesBC);

% Initialize variables
modelTextBCStructure = [];
% Get IQMstructure
IQMstructure = IQMstruct(iqm);
% Parse structure into the modelTextBCStructure description

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define filename used for intermediate saving
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = strcat(tempdirIQM,'tempsavingfile.temp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelTextBCStructure.name = IQMstructure.name;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelTextBCStructure.notes = IQMstructure.notes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS THE states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(filename,'w');
% first the definition of the states+ODEs that need to be described by ODEs
for k = 1:length(stateNamesODE),
    % get index of state name
    index = strmatchIQM(stateNamesODE{k},{IQMstructure.states.name},'exact');
    fprintf(fid,'d/dt(%s) = %s\n',stateNamesODE{k},IQMstructure.states(index).ODE);
end
if ~isempty(stateNamesODE),
	fprintf(fid,'\n');
end

% WRITE OUT ALGEBRAIC RULES
for k = 1:length(IQMstructure.algebraic),
    if ~isempty(IQMstructure.algebraic(k).name),
        fprintf(fid,'0 = %s : %s\n',IQMstructure.algebraic(k).formula,IQMstructure.algebraic(k).name);
    else
        fprintf(fid,'0 = %s\n',IQMstructure.algebraic(k).formula);
    end
end
if ~isempty(IQMstructure.algebraic),
	fprintf(fid,'\n');
end

% WRITE OUT INITIAL CONDITIONS
% definition of the initial conditions of the states.
% additionally for each component optional additional information
% is processed and written behind the initial conditions.
% states.type
% states.compartment
% states.unittype
% now construct the states text
for k = 1:length(IQMstructure.states),
    name = IQMstructure.states(k).name;
    if isnumeric(IQMstructure.states(k).initialCondition),
        fprintf(fid,'%s(0) = %1.6g',name,IQMstructure.states(k).initialCondition);
    else
        fprintf(fid,'%s(0) = %s',name,IQMstructure.states(k).initialCondition);
    end
    % construct and the additional information
    type = strtrim(IQMstructure.states(k).type);
    compartment = strtrim(IQMstructure.states(k).compartment);
    unittype = strtrim(IQMstructure.states(k).unittype);
    % check if all information is given
    if ~isempty(strfind(lower(type),'specie')),
        informationText = strcat('{',type,':',compartment,':',unittype,'}');
    elseif ~isempty(strfind(lower(type),'compartment')),
        informationText = strcat('{',type,':',compartment,'}');
    elseif ~isempty(strfind(lower(type),'parameter')),
        informationText = strcat('{',type,'}');
    else
        informationText = '';
    end
    fprintf(fid,' %s',informationText);
    if ~isempty(IQMstructure.states(k).notes)
        fprintf(fid,' %% %s',IQMstructure.states(k).notes);
    end
    fprintf(fid,'\n');
end
% do the same for algebraic variable initial conditions
for k = 1:length(IQMstructure.algebraic),
    if ~isempty(IQMstructure.algebraic(k).name),
        name = IQMstructure.algebraic(k).name;
        if isnumeric(IQMstructure.algebraic(k).initialCondition),
            fprintf(fid,'%s(0) = %1.6g',name,IQMstructure.algebraic(k).initialCondition);
        else
            fprintf(fid,'%s(0) = %s',name,IQMstructure.algebraic(k).initialCondition);
        end
        % construct the additional information
        type = strtrim(IQMstructure.algebraic(k).type);
        compartment = strtrim(IQMstructure.algebraic(k).compartment);
        unittype = strtrim(IQMstructure.algebraic(k).unittype);
        % check if all information is given
        if ~isempty(strfind(lower(type),'specie')),
            informationText = strcat('{',type,':',compartment,':',unittype,'}');
        elseif ~isempty(strfind(lower(type),'compartment')),
            informationText = strcat('{',type,':',compartment,'}');
        elseif ~isempty(strfind(lower(type),'parameter')),
            informationText = strcat('{',type,'}');
        else
            informationText = '';
        end
        fprintf(fid,' %s',informationText);
        if ~isempty(IQMstructure.algebraic(k).notes)
            fprintf(fid,' %% %s',IQMstructure.algebraic(k).notes);
        end
        fprintf(fid,'\n');
    end
end
% rewind the file and read it out
fclose(fid); fid = fopen(filename,'r');
modelTextBCStructure.states = fread(fid,inf,'uint8=>char')';
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add the parameters
% additionally for each parameter optional additional information
% is processed and written behind the parameters.
% parameters.type
% parameters.compartment
% parameters.unittype
fid = fopen(filename,'w');
for k = 1:length(IQMstructure.parameters),
    fprintf(fid,'%s = %1.6g',IQMstructure.parameters(k).name,IQMstructure.parameters(k).value);
    % construct and the additional information
    type = strtrim(IQMstructure.parameters(k).type);
    compartment = strtrim(IQMstructure.parameters(k).compartment);
    unittype = strtrim(IQMstructure.parameters(k).unittype);
    % check if all information is given
    if ~isempty(strfind(lower(type),'specie')),
        informationText = strcat('{',type,':',compartment,':',unittype,'}');
    elseif ~isempty(strfind(lower(type),'compartment')),
        informationText = strcat('{',type,':',compartment,'}');
    elseif ~isempty(strfind(lower(type),'parameter')),
        informationText = strcat('{',type,'}');
    else
        informationText = '';
    end
    fprintf(fid,' %s',informationText);
    if ~isempty(IQMstructure.parameters(k).notes)
        fprintf(fid,' %% %s',IQMstructure.parameters(k).notes);
    end
    fprintf(fid,'\n');
end
% rewind the file and read it out
fclose(fid); fid = fopen(filename,'r');
modelTextBCStructure.parameters = fread(fid,inf,'uint8=>char')';
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(filename,'w');
for k = 1:length(IQMstructure.variables),
    fprintf(fid,'%s = %s',IQMstructure.variables(k).name,IQMstructure.variables(k).formula);
    % construct and the additional information
    type = strtrim(IQMstructure.variables(k).type);
    compartment = strtrim(IQMstructure.variables(k).compartment);
    unittype = strtrim(IQMstructure.variables(k).unittype);
    % check if all information is given
    if ~isempty(strfind(lower(type),'specie')),
        informationText = strcat('{',type,':',compartment,':',unittype,'}');
    elseif ~isempty(strfind(lower(type),'compartment')),
        informationText = strcat('{',type,':',compartment,'}');
    elseif ~isempty(strfind(lower(type),'parameter')),
        informationText = strcat('{',type,'}');
    else
        informationText = '';
    end
    fprintf(fid,' %s',informationText);
    if ~isempty(IQMstructure.variables(k).notes)
        fprintf(fid,' %% %s',IQMstructure.variables(k).notes);
    end
    fprintf(fid,'\n');
end
% rewind the file and read it out
fclose(fid); fid = fopen(filename,'r');
modelTextBCStructure.variables = fread(fid,inf,'uint8=>char')';
fclose(fid);

[N, stateNames] = reverseMoietyConservationsIQM(iqm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS THE REACTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use the stoichiometric matrix determined above (by eventually adding
% species that are defined by variables)
% get reaction names, kinetics, and reversibility flag
[reactionNames,reactionFormulas,reactionRevFlags] = IQMreactions(iqm);
% now the reaction text can be built
fid = fopen(filename,'w');
% cycle through the columns of the stoichiometric matrix N to build the 
% reaction expressions
for k1 = 1:size(N,2),
    Ncol = N(:,k1);
    % first get the substrates by finding the negative elements in Ncol
    substrateIndices = find(Ncol < 0);
    % then get the products by finding the positive elements in Ncol
    productIndices = find(Ncol > 0);
    % determine the needed information
    reactionName = reactionNames{k1};
    reactionRevFlag = reactionRevFlags(k1);
    reactionFormula = reactionFormulas{k1};
    substrateNames = stateNames(substrateIndices);
    productNames = stateNames(productIndices);
    substrateStoichiometries = abs(Ncol(substrateIndices));
    productStoichiometries = abs(Ncol(productIndices));
    % if reversible split up the reaction rate in two parts. if this is not
    % possible, issue a warning and set the reaction as irreversible
    if reactionRevFlag ~= 0,
        irreversibleRates = explodePCIQM(reactionFormula,'-');
        if length(irreversibleRates) ~= 2,
            % Need to have two parts that are separated by a '-' sign. (Forward
            % first, then reverse reaction kinetics).
            warning(sprintf('Reaction ''%s'' is marked reversible but rate expression is not given in the\nrequired format (R = Rf-Rr). The reaction will be treated as irreversible.\n', reactionName));
            reactionRevFlag = 0;
        else
            reactionForward = irreversibleRates{1};
            reactionReverse = irreversibleRates{2};
        end
    end
    % format the output of the reaction text
    % first the reaction expression, e.g. 2*A + 4*C => 3*B
    % the substrates
    if ~isempty(substrateNames),
        if substrateStoichiometries(1) ~= 1,
            fprintf(fid,'%g*%s',substrateStoichiometries(1),substrateNames{1});
        else
            fprintf(fid,'%s',substrateNames{1});
        end
    end
    for k2 = 2:length(substrateNames),
        if substrateStoichiometries(k2) ~= 1,
            fprintf(fid,'+%g*%s',substrateStoichiometries(k2),substrateNames{k2});
        else
            fprintf(fid,'+%s',substrateNames{k2});
        end
    end
    % the reaction equation sign
    if reactionRevFlag == 0,
        fprintf(fid,' => ');  % irreversible
    else
        fprintf(fid,' <=> '); % reversible
    end
    % the products
    if ~isempty(productNames),
        if productStoichiometries(1) ~= 1,
            fprintf(fid,'%g*%s',productStoichiometries(1),productNames{1});
        else
            fprintf(fid,'%s',productNames{1});
        end
    end
    for k2 = 2:length(productNames),
        if productStoichiometries(k2) ~= 1,
            fprintf(fid,'+%g*%s',productStoichiometries(k2),productNames{k2});
        else
            fprintf(fid,'+%s',productNames{k2});
        end
    end
    % separator and reaction name
    fprintf(fid,' : %s',reactionName); 
    % fast flag
    if IQMstructure.reactions(k1).fast == 1,
        fprintf(fid,' {fast}'); 
    end
    % notes
    if ~isempty(IQMstructure.reactions(k1).notes),
        fprintf(fid,' %% %s',IQMstructure.reactions(k1).notes);
    end
    % new line 
    fprintf(fid,'\n'); 
    % now the reaction rate expression(s)
    if reactionRevFlag == 0,
        fprintf(fid,'\tvf = %s\n',reactionFormula); 
    else
        fprintf(fid,'\tvf = %s\n\tvr = %s\n',reactionForward,reactionReverse); 
    end
    % new line after reaction definition
    fprintf(fid,'\n'); 
end
% rewind the file and read it out
fclose(fid); fid = fopen(filename,'r');
modelTextBCStructure.reactions = fread(fid,inf,'uint8=>char')';
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(filename,'w');
for k = 1:length(IQMstructure.functions),
    fprintf(fid,'%s(%s) = %s',IQMstructure.functions(k).name,IQMstructure.functions(k).arguments,IQMstructure.functions(k).formula);
    if ~isempty(IQMstructure.functions(k).notes)
        fprintf(fid,' %% %s\n',IQMstructure.functions(k).notes);
    else
        fprintf(fid,'\n');
    end
end
% rewind the file and read it out
fclose(fid); fid = fopen(filename,'r');
modelTextBCStructure.functions = fread(fid,inf,'uint8=>char')';
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Events
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(filename,'w');
for k = 1:length(IQMstructure.events),
    fprintf(fid,'%s = %s',IQMstructure.events(k).name,IQMstructure.events(k).trigger);
    for k2 = 1:length(IQMstructure.events(k).assignment),
        fprintf(fid,',%s,%s',IQMstructure.events(k).assignment(k2).variable,IQMstructure.events(k).assignment(k2).formula);
    end
    if ~isempty(IQMstructure.events(k).notes)
        fprintf(fid,' %% %s\n',IQMstructure.events(k).notes);
    else
        fprintf(fid,'\n');
    end
end
% rewind the file and read it out
fclose(fid); fid = fopen(filename,'r');
modelTextBCStructure.events = fread(fid,inf,'uint8=>char')';
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelTextBCStructure.functionsMATLAB = IQMstructure.functionsMATLAB;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DELETE WHITESPACES IN STRINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Useful for taking away whitespaces in kineticLaw formulas, as
% seen in some example models
function [outputString] = removeWhiteSpace(inputString)
outputString = strrep(inputString,' ','');
return
