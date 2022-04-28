function [IQMstructure,errorMsg] = convertSBML2IQM(SBMLmodel,SBMLmodelFile)
% convertSBML2IQM
% Converting a SBML Level 2 MATLAB structure to the object model structure
% used in the toolbox
%
% Private method 
%
% [IQMstructure,errorMsg] = convertSBML2IQM(SBMLmodel);
%
% Simple error checking is provided. In the case of an error an empty
% structure is returned.

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE ERROR MESSAGES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
errorMsg = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do some checks if supported or not
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
contents = fileread(SBMLmodelFile);

% Stoichiometry math wrongly handeled in libSBML
if ~isempty(strfind(contents,'<stoichiometryMath>')),
    errorMsg = sprintf('%s\nStoichiometry math construct is wrongly handeled in libSBML and thus not supported currently in IQM Tools.\n',errorMsg);
end

% Do not consider delay things in the model
if ~isempty(strfind(contents,'<delay>')),
    errorMsg = sprintf('%s\nDelays in SBML models are not supported in IQM Tools.\n',errorMsg);
end

% Do not consider delay things in the model
if ~isempty(strfind(contents,'symbols/delay')),
    errorMsg = sprintf('%s\nDelays in SBML models are not supported in IQM Tools.\n',errorMsg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start conversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global SBMLtimesymbol

% additional global variable to record species that have amount initial
% units but are converted to concentration ... and the corresponding
% compartment
global amount2concentration
amount2concentration = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE AN EMPTY IQMSTRUCTURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IQMstructure = struct(IQMmodel());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK TIME SYMBOL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From SBML Toolbobx 2.0.2 on a time_symbol field exists in the
% SBMLmodel. We need to get the time symbol if it is there and replace
% all occurrences in the model with 'time'.
if isfield(SBMLmodel,'time_symbol'),
    SBMLtimesymbol = SBMLmodel.time_symbol;
else
    SBMLtimesymbol = '';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NAME AND NOTES (INCLUDE IN STRUCTURE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp(SBMLmodel.name,''),
    IQMstructure.name = SBMLmodel.name;
elseif ~strcmp(SBMLmodel.id,''),
    IQMstructure.name = SBMLmodel.id;
else
    IQMstructure.name = 'Unknown model';
end
% Indicate in notes that the model has been created by import from SBML.
IQMstructure.notes = strtrim(convert2IQMNotes(SBMLmodel.notes,0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS (INCLUDE IN STRUCTURE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:length(SBMLmodel.functionDefinition),
    functionName = SBMLmodel.functionDefinition(k).id;
    functionMath = SBMLmodel.functionDefinition(k).math;
    % take away the lambda( ... ) and only leave the ...
    expression = removeWhiteSpace(functionMath(8:length(functionMath)-1));
    % take away power functions
    expression = exchangepowerexp(expression);
    % explode at the comma (not within parentheses)
    [elements] = explodePCIQM(expression);
    % last element is the formula, the others are the arguments
    functionFormula = elements{end};
    functionArguments = expression(1:end-length(functionFormula)-1);
    % replace eventual MathML expressions by MATLAB expressions
    [functionFormula] = replaceMathMLexpressions(functionFormula);
    % include function into structure
    
    IQMstructure.functions(k).name = functionName;
    IQMstructure.functions(k).arguments = functionArguments;
    IQMstructure.functions(k).formula = functionFormula;
    IQMstructure.functions(k).notes = convert2IQMNotes(SBMLmodel.functionDefinition(k).notes,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOME COUNTER DEFINITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% States and parameters are included in the structure in their order of
% appearance. For variables no counters are needed, since their order 
% is defined by the ordering of the scalar rules and not by their order of 
% appearance (see below)
numberStates = 1;
numberParameters = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ORDERING OF VARIABLES (SCALAR RULES)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scalar rules lead always to the definition of variables. The order of
% appearance of the rules is important and thus the order of appearance of
% the variables should reflect the ordering of the rules. We determine an
% array in which the the index i determines the priority of the scalar rule, and
% the i-th array element the index of the corresponding scalar rule. Rate
% rules do not have to be considered since they are evaluated at last
% anyway
orderScalarRules = strmatchIQM('SBML_ASSIGNMENT_RULE',{SBMLmodel.rule.typecode});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECIES (INCLUDE IN STRUCTURE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Species of the model are included in the structure as follows:
%
% Non-boundary species:
%   - rate rule: as state
%   - scalar rule: as variable
%   - no rule: as state
%
% Boundary species:
%   - rate rule: as states
%   - scalar rule: as variable
%   - no rule: as parameter
%
% If initial values of species given as concentrations:
%   - They will be treated as concentrations, this means the ODE will be
%     adjusted by diving the amount rate with the corresponding compartment
%     volume.
%
% If initial values of species given as amounts:
%   - They will be treated as amounts, this means the ODE will NOT be
%     adjusted!
%
% The rate formulas in the kinetic laws are just taken as they are. This
% means a mix of concentration and amount species might occur - the user 
% has him/herself to take care to have the correct rate laws!
%
% Species that have amount units will be divided with their compartment volume in each expression later on
speciesAmount_SAVEnames = {};
speciesAmount_SAVEcompartments = {};
%
% Species that are defined by rate rules need this rate rule to be multiplied by the compartment the species is in.

%
speciesODElist = [];
speciesODElistIndex = 1;
% Check if there are any species. 
if isempty(SBMLmodel.species),
%    errorMsg = sprintf('%s\n%s\n',errorMsg,'No species are defined in the SBML model');
%     warning('No species are defined in the SBML model - uncertain if model conversion will lead to a working IQMmodel!');
end
% Cycle through all the species and determine how they are to be included
% in the model structure
for k1 = 1:length(SBMLmodel.species),
%%
    speciesName = SBMLmodel.species(k1).id;
    speciesCompartmentName = SBMLmodel.species(k1).compartment;
    speciesNotes = convert2IQMNotes(SBMLmodel.species(k1).notes,1);
    % get the compartment size for this species
    makeItAmountAnyway = 0;
    for k2 = 1:length(SBMLmodel.compartment),
        if strcmp(speciesCompartmentName, SBMLmodel.compartment(k2).id),
            compartmentSize = SBMLmodel.compartment(k2).size;
            if isempty(compartmentSize),
                errorMsg = sprintf('%s\nNo compartment size defined for compartment ''%s''\n',errorMsg,speciesCompartmentName);
                compartmentSize = 1; % set to whatever, error will still occurr!
            end
            spatialDimensions = SBMLmodel.compartment(k2).spatialDimensions;
            if spatialDimensions == 0 || compartmentSize == 0,
                makeItAmountAnyway = 1;
            end
            break;
        end
    end
    % determine the initial value for species
    % if hasOnlySubstanceUnits is true then it is amount, otherwise it is concentration
    initialCondition = [];
    if SBMLmodel.species(k1).hasOnlySubstanceUnits || makeItAmountAnyway || (SBMLmodel.species(k1).isSetInitialAmount && ~SBMLmodel.species(k1).hasOnlySubstanceUnits ),% && ~isempty(SBMLmodel.species(k1).substanceUnits)),
        % It is amount ... need to calculate if needed
        % Only if substanceunits is not substance
%         if ~strcmp(SBMLmodel.species(k1).substanceUnits,'substance') && ~SBMLmodel.species(k1).hasOnlySubstanceUnits,
            speciesAmount_SAVEnames{end+1} = SBMLmodel.species(k1).id;
            speciesAmount_SAVEcompartments{end+1} = SBMLmodel.species(k1).compartment;
%         end

        speciesUnits = 'amount';
        if SBMLmodel.species(k1).isSetInitialAmount,
            initialCondition = SBMLmodel.species(k1).initialAmount;
        elseif SBMLmodel.species(k1).isSetInitialConcentration,
            initialCondition = SBMLmodel.species(k1).initialConcentration*compartmentSize;
        else 
            initialCondition = NaN;
        end
    else
        % It is concentration
        speciesUnits = 'concentration';
        if SBMLmodel.species(k1).isSetInitialAmount,
            initialCondition = SBMLmodel.species(k1).initialAmount/compartmentSize;
        elseif SBMLmodel.species(k1).isSetInitialConcentration,
            initialCondition = SBMLmodel.species(k1).initialConcentration;
        else 
            initialCondition = NaN;
        end
    end
% initialCondition
% speciesUnits
    %%
    % Check the number of rules the species is in and get the index
    % of the last rule the species is used in together with its type
    % ('rate' or 'scalar') and the right hand side expression 'formula'
    [numberRulesSpecies,lastRuleIndex,lastRuleType,lastRuleFormula,errorMsg] = getRules(speciesName,SBMLmodel,errorMsg);
    if numberRulesSpecies > 1,
        % Species defined by more than one rule => error
        errorMsg = sprintf('%s\nSpecies ''%s'' defined by more than one rule\n',errorMsg,speciesName);
    end
    % Process all the species that DO NOT have a boundary condition
    % In case the species is non-constant it is included as for Level 1
    % In case the species is constant it is included as ODE with RHS 0 and initial condition as value
    if SBMLmodel.species(k1).constant == 0,
        if SBMLmodel.species(k1).boundaryCondition == 0,
            if numberRulesSpecies == 0,
                % add species to the ODE list (for use below in ODE
                % construction)
                speciesODElist(speciesODElistIndex).name = speciesName;
                speciesODElist(speciesODElistIndex).stateIndex = numberStates;
                speciesODElist(speciesODElistIndex).units = speciesUnits;
                speciesODElistIndex = speciesODElistIndex + 1;
                % Include species as state
                IQMstructure.states(numberStates).name = speciesName;
                % Check if an initial value is set for a species
                if isempty(initialCondition),
%                     errorMsg = sprintf('%s\nNo initial value defined for species ''%s''\n',errorMsg,speciesName);
                    initialCondition = 0;
                end
                IQMstructure.states(numberStates).initialCondition = initialCondition;
                % Note that this state is a 'species'
                % Notes are not used by the toolbox for other things as
                % documentation
                IQMstructure.states(numberStates).notes = speciesNotes;
                IQMstructure.states(numberStates).type = 'isSpecie';
                IQMstructure.states(numberStates).compartment = SBMLmodel.species(k1).compartment;
                IQMstructure.states(numberStates).unittype = speciesUnits;                
                % Initialize the ODE field for this state
                IQMstructure.states(numberStates).ODE = '';
                % Increment
                numberStates = numberStates+1;
            elseif numberRulesSpecies == 1,
                % Before including the species in the model structure we have
                % to check if this species is also used in reactions as product
                % or reactant. In this case an error message will be issued,
                % since non-boundary species are not allowed to appear in
                % reactions and rules.
                reactionSpeciesPresent = checkReactionsForSpecies(speciesName,SBMLmodel);
                if reactionSpeciesPresent,
                    % Species is also present in a reaction. This means the SBML
                    % model is inconsistent. => error
                    errorMsg = sprintf('%s\nSpecies ''%s'' defined by rule and altered by reactions\n',errorMsg,speciesName);
                end
                % Even in the case of an error we continue here with the
                % processing - The error message will take care of returning
                % an empty model structure at the end
                if strcmp(lastRuleType,'rate'),
                    % Species added as a state, since rule of type 'rate'
                    IQMstructure.states(numberStates).name = speciesName;
                    % check if an initial value is set for a species
                    if isempty(initialCondition),
%                         errorMsg = sprintf('%s\nNo initial value defined for species ''%s''\n',errorMsg,speciesName);
                        initialCondition = 0;
                    end
                    % Include initial value into structure
                    IQMstructure.states(numberStates).initialCondition = initialCondition;
                    % Include the RHS formula as the ODE for this species state
                    % in case of speciesUnits='concentration' keep as is the rule.
                    % But in case of amount then multiply by compartment.
                    
                    formula_use = lastRuleFormula;
                    if strcmp(speciesUnits,'amount'),
                        formula_use = ['(' lastRuleFormula ')*' SBMLmodel.species(k1).compartment];
                    end
                    IQMstructure.states(numberStates).ODE = formula_use;
                    
                    % Set a note that type of state is 'species'
                    IQMstructure.states(numberStates).notes = speciesNotes;
                    % initialize type, compartment, and unittype fields
                    IQMstructure.states(numberStates).type = 'isSpecie';
                    IQMstructure.states(numberStates).compartment = SBMLmodel.species(k1).compartment;
                    IQMstructure.states(numberStates).unittype = speciesUnits;                
                    % Increment
                    numberStates = numberStates+1;
                else
                    % Species added as a variable, since rule of type 'scalar'
                    % Determine the index where to include the variable (based
                    % on the ordering of the scalar rules)
                    
                    formula_use = lastRuleFormula;
                    if sum([SBMLmodel.species.isSetInitialConcentration]) == 0,
                        formula_use = ['(' lastRuleFormula ')*' SBMLmodel.species(k1).compartment];
                    end
                    
                    indexVariable = find(orderScalarRules==lastRuleIndex);
                    IQMstructure.variables(indexVariable).name = speciesName;
                    IQMstructure.variables(indexVariable).formula = formula_use;
                    IQMstructure.variables(indexVariable).notes = speciesNotes;
                    % initialize type, compartment, and unittype fields
                    IQMstructure.variables(indexVariable).type = 'isSpecie';
                    IQMstructure.variables(indexVariable).compartment = SBMLmodel.species(k1).compartment;
                    IQMstructure.variables(indexVariable).unittype = speciesUnits;
                end
            end
            % The case where numberRulesSpecies > 1 does not have to be treated
            % here since already an error message is set above. This error
            % message will lead to an empty model structure at the end. Even
            % if the error is set the rest of the model will be processed.
            % This has the advantage that all errors can be
            % detected at once and makes the control structure of this script
            % simpler.
        else
            % Process all the species that DO have a boundary condition
            if numberRulesSpecies == 0,
                % Include boundary species as a parameter
                IQMstructure.parameters(numberParameters).name = speciesName;
                % Check if an initial amount is set for a species
                if isempty(initialCondition),
%                     errorMsg = sprintf('%s\nNo initial value defined for species ''%s''\n',errorMsg,speciesName);
                    initialCondition = 0;
                end
                % Include initial value
                IQMstructure.parameters(numberParameters).value = initialCondition;
                % Set note that type of the parameter is 'boundaey species'
                IQMstructure.parameters(numberParameters).notes = speciesNotes;
                % initialize type, compartment, and unittype fields
                IQMstructure.parameters(numberParameters).type = 'isSpecie';
                IQMstructure.parameters(numberParameters).compartment = SBMLmodel.species(k1).compartment;
                IQMstructure.parameters(numberParameters).unittype = speciesUnits;
                % Increment
                numberParameters = numberParameters+1;
            elseif numberRulesSpecies == 1,
                % Boundary species may appear as products and reactants in
                % reactions.
                if strcmp(lastRuleType,'rate'),
                    % Boundary species added as a state, since rule of type 'rate'
                    IQMstructure.states(numberStates).name = speciesName;
                    % check if an initial amount is set for a species
                    if isempty(initialCondition),
%                         errorMsg = sprintf('%s\nNo initial value defined for species ''%s''\n',errorMsg,speciesName);
                        initialCondition = 0;
                    end
                    % Include initial value in structure
                    IQMstructure.states(numberStates).initialCondition = initialCondition;
                    % Include the RHS formula as the ODE for this species state
                    IQMstructure.states(numberStates).ODE = lastRuleFormula;
                    % Set note that type of the state is 'boundary species'
                    IQMstructure.states(numberStates).notes = speciesNotes;
                    % initialize type, compartment, and unittype fields
                    IQMstructure.states(numberStates).type = 'isSpecie';
                    IQMstructure.states(numberStates).compartment = SBMLmodel.species(k1).compartment;
                    IQMstructure.states(numberStates).unittype = speciesUnits;
                    % Increment
                    numberStates = numberStates+1;
                else
                    % Boundary species added as a variable, since rule of type 'scalar'
                    % Determine the index where to include the variable (based
                    % on the ordering of the scalar rules)
                    indexVariable = find(orderScalarRules==lastRuleIndex);
                    IQMstructure.variables(indexVariable).name = speciesName;
                    IQMstructure.variables(indexVariable).formula = lastRuleFormula;
                    IQMstructure.variables(indexVariable).notes = speciesNotes;
                    % initialize type, compartment, and unittype fields
                    IQMstructure.variables(indexVariable).type = 'isSpecie';
                    IQMstructure.variables(indexVariable).compartment = SBMLmodel.species(k1).compartment;
                    IQMstructure.variables(indexVariable).unittype = speciesUnits;
                end
            end
        end
    else
        % Species defined as being constant!
        % Thus include it as state with ODE RHS 0
        IQMstructure.states(numberStates).name = speciesName;
        % check if an initial amount is set for a species
        if isempty(initialCondition),
            %                         errorMsg = sprintf('%s\nNo initial value defined for species ''%s''\n',errorMsg,speciesName);
            initialCondition = 0;
        end
        % Include initial value in structure
        IQMstructure.states(numberStates).initialCondition = initialCondition;
        % 0 RHS since constant
        IQMstructure.states(numberStates).ODE = '0';
        % Set note that type of the state is 'boundary species'
        IQMstructure.states(numberStates).notes = speciesNotes;
        % initialize type, compartment, and unittype fields
        IQMstructure.states(numberStates).type = 'isSpecie';
        IQMstructure.states(numberStates).compartment = SBMLmodel.species(k1).compartment;
        IQMstructure.states(numberStates).unittype = speciesUnits;
        % Increment
        numberStates = numberStates+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS (INCLUDE IN STRUCTURE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First process the global parameters and check if rules exist for them
%   - rate rule: as state
%   - scalar rule: as variable
%   - no rule: as parameter
%
% Then process the local parameters. They can only be parameters, no 
% rules can exist for them. In case of collisions:
%
% parameter-parameter: colliding parameters are prefixed by their reaction 
%                      names (global parameters are not prefixed)
%
% parameter-species: are prefixed by reaction name
%
% parameter-compartment: are prefixed by reactio name
%
% Start with global parameters
for k = 1:length(SBMLmodel.parameter),
    parameterName = SBMLmodel.parameter(k).id;
    parameterValue = SBMLmodel.parameter(k).value;
    parameterNotes = convert2IQMNotes(SBMLmodel.parameter(k).notes,1);
    % Check the number of rules the parameter is in and get the index
    % of the last rule the parameter is used in together with its type
    % ('rate' or 'scalar') and the right hand side expression 'formula'
    [numberRulesParameter,lastRuleIndex,lastRuleType,lastRuleFormula,errorMsg] = getRules(parameterName,SBMLmodel,errorMsg);
    if numberRulesParameter > 1,
        % Parameter defined by more than one rule => error
        errorMsg = sprintf('%s\nParameter ''%s'' defined by more than one rule\n',errorMsg,parameterName);
    end
    if numberRulesParameter == 0 || SBMLmodel.parameter(k).constant,
        % Include parameter as parameter
        IQMstructure.parameters(numberParameters).name = parameterName;
        % Check if a value is set for the parameter
        if length(parameterValue)==0,
            errorMsg = sprintf('%s\nNo value defined defined for parameter ''%s''\n',errorMsg,parameterName);
        end
        % Include value in IQMstructure
        IQMstructure.parameters(numberParameters).value = parameterValue;
        % Set note that type of the parameter is 'global'
        IQMstructure.parameters(numberParameters).notes = parameterNotes; 
        % initialize type, compartment, and unittype fields
        IQMstructure.parameters(numberParameters).type = 'isParameter';
        IQMstructure.parameters(numberParameters).compartment = '';
        IQMstructure.parameters(numberParameters).unittype = '';                
        % Increment
        numberParameters = numberParameters+1;
    elseif numberRulesParameter == 1,
        % One rule has been detected for this parameter
        if strcmp(lastRuleType,'rate'),
            % Parameter added as a state, since rule of type 'rate'
            IQMstructure.states(numberStates).name = parameterName;
            % check if an initial amount is set for the parameter
            if length(parameterValue)==0,
                errorMsg = sprintf('%s\nNo value defined for parameter ''%s''\n',errorMsg,parameterName);
            end
            % Include value as initial value in IQMstructure
            IQMstructure.states(numberStates).initialCondition = parameterValue;
            % Include the RHS formula as the ODE for this parameter state
            IQMstructure.states(numberStates).ODE = lastRuleFormula;
            % Set note that type of the state is 'parameter rule'
            IQMstructure.states(numberStates).notes = parameterNotes;
            % initialize type, compartment, and unittype fields
            IQMstructure.states(numberStates).type = 'isParameter';
            IQMstructure.states(numberStates).compartment = '';
            IQMstructure.states(numberStates).unittype = '';                
            % Increment
            numberStates = numberStates+1;
        else
            % Parameter added as a variable, since rule of type 'scalar'
            % Determine the index where to include the variable (based
            % on the ordering of the scalar rules)
            indexVariable = find(orderScalarRules==lastRuleIndex);
            IQMstructure.variables(indexVariable).name = parameterName;
            IQMstructure.variables(indexVariable).formula = lastRuleFormula;
            IQMstructure.variables(indexVariable).notes = parameterNotes;
            % initialize type, compartment, and unittype fields
            IQMstructure.variables(indexVariable).type = 'isParameter';
            IQMstructure.variables(indexVariable).compartment = '';
            IQMstructure.variables(indexVariable).unittype = '';               
        end
    end
end
% Then include the local parameters within the different reactions
% Name collisions are avoided by adding the reaction name in front of each
% local parameter.
for k1 = 1:length(SBMLmodel.reaction),
    % get the current reaction name
    reactionName = SBMLmodel.reaction(k1).id;
    if ~isempty(SBMLmodel.reaction(k1).kineticLaw),
        for k2 = 1:length(SBMLmodel.reaction(k1).kineticLaw.parameter),
            % get the current local parameter in reaction reactionName
            parameterName = SBMLmodel.reaction(k1).kineticLaw.parameter(k2).id;
            parameterValue = SBMLmodel.reaction(k1).kineticLaw.parameter(k2).value;
            parameterNotes = convert2IQMNotes(SBMLmodel.reaction(k1).kineticLaw.parameter(k2).notes,1);
            % check if a value is set for the parameter
            if length(parameterValue)==0,
                errorMsg = sprintf('%s\nNo value defined for parameter ''%s'' in reaction ''%s''\n',errorMsg,parameterName,reactionName);
            end
            % add the reaction name to the parameter and update the kinetic rate law
            parameterName = char([double(reactionName) double('_') double(parameterName)]);  % fastest strcat ;)
            % Include parameter in the structure
            IQMstructure.parameters(numberParameters).name = parameterName;
            IQMstructure.parameters(numberParameters).value = parameterValue;
            IQMstructure.parameters(numberParameters).notes = parameterNotes;
            % initialize type, compartment, and unittype fields
            IQMstructure.parameters(numberParameters).type = 'isParameter';
            IQMstructure.parameters(numberParameters).compartment = '';
            IQMstructure.parameters(numberParameters).unittype = '';
            numberParameters = numberParameters+1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARTMENTS (INCLUDE IN STRUCTURE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Include compartments
% No rule => as parameters
%
% Scalar rule => as variables
%       
% Rate rule => as states
%
for k = 1:length(SBMLmodel.compartment),
    compartmentName = SBMLmodel.compartment(k).id;
    compartmentValue = SBMLmodel.compartment(k).size;
    compartmentNotes = convert2IQMNotes(SBMLmodel.compartment(k).notes,1);
    % Check the number of rules the compartment is in and get the index
    % of the last rule the compartment is used in together with its type
    % ('rate' or 'scalar') and the right hand side expression 'formula'
    [numberRulesCompartment,lastRuleIndex,lastRuleType,lastRuleFormula,errorMsg] = getRules(compartmentName,SBMLmodel,errorMsg);
    if numberRulesCompartment > 1,
        % Compartment defined by more than one rule => error
        errorMsg = sprintf('%s\nCompartment ''%s'' defined by more than one rule\n',errorMsg,compartmentName);
    end
    if numberRulesCompartment == 1 && SBMLmodel.compartment(k).constant,
        % compartment defined as constant, but rule defined for it
        errorMsg = sprintf('%s\nCompartment ''%s'' defined as having constant size but rule exists for it.\n',errorMsg,compartmentName);
    end
    if numberRulesCompartment == 0 || SBMLmodel.compartment(k).constant,
        % check compartment size and return an error in case it is "NaN"
        if isnan(compartmentValue),
%            errorMsg = sprintf('%s\nNo size for compartment ''%s'' given\n',errorMsg,compartmentName);
             compartmentValue = 1;   
        end
        % Include compartment as parameter
        IQMstructure.parameters(numberParameters).name = compartmentName;
        % Check if a value is set for the parameter
        if isempty(compartmentValue),
%             errorMsg = sprintf('%s\nNo value defined defined for compartment ''%s''\n',errorMsg,compartmentName);
            compartmentValue = 1;
        end
        % Include value in IQMstructure
        IQMstructure.parameters(numberParameters).value = compartmentValue;
        % Set note that type of the parameter is 'compartment size'
        IQMstructure.parameters(numberParameters).notes = compartmentNotes;
        % initialize type, compartment, and unittype fields
        IQMstructure.parameters(numberParameters).type = 'isCompartment';
        IQMstructure.parameters(numberParameters).compartment = SBMLmodel.compartment(k).outside;
        IQMstructure.parameters(numberParameters).unittype = '';
        % Increment
        numberParameters = numberParameters+1;
    elseif numberRulesCompartment == 1,
        % One rule has been detected for this compartment
        if strcmp(lastRuleType,'rate'),
            % check compartment size and return an error in case it is "NaN"
            if isnan(compartmentValue),
%                 errorMsg = sprintf('%s\nNo size for compartment ''%s'' given\n',errorMsg,compartmentName);
                compartmentValue = 1;
            end
            % Compartment added as a state, since rule of type 'rate'
            IQMstructure.states(numberStates).name = compartmentName;
            % check if an initial amount is set for the compartment
            if isempty(compartmentValue),
%                 errorMsg = sprintf('%s\nNo value defined for compartment ''%s''\n',errorMsg,compartmentName);
                compartmentValue = 1;
            end
            % Include value as initial value in IQMstructure
            IQMstructure.states(numberStates).initialCondition = compartmentValue;
            % Include the RHS formula as the ODE for this compartment state
            IQMstructure.states(numberStates).ODE = lastRuleFormula;
            % Set note that type of the state is 'compartment size'
            IQMstructure.states(numberStates).notes = compartmentNotes;
            % initialize type, compartment, and unittype fields
            IQMstructure.states(numberStates).type = 'isCompartment';
            IQMstructure.states(numberStates).compartment = SBMLmodel.compartment(k).outside;
            IQMstructure.states(numberStates).unittype = '';
            % Increment
            numberStates = numberStates+1;
            % compartment_rate variable right hand side
            compartmentDotRHS = lastRuleFormula;
        else
            % Compartment added as a variable, since rule of type 'scalar'
            % Determine the index where to include the variable (based
            % on the ordering of the scalar rules)
            indexVariable = find(orderScalarRules==lastRuleIndex);
            IQMstructure.variables(indexVariable).name = compartmentName;
            IQMstructure.variables(indexVariable).formula = lastRuleFormula;
            IQMstructure.variables(indexVariable).notes = compartmentNotes;
            % initialize type, compartment, and unittype fields
            IQMstructure.variables(indexVariable).type = 'isCompartment';
            IQMstructure.variables(indexVariable).compartment = SBMLmodel.compartment(k).outside;
            IQMstructure.variables(indexVariable).unittype = '';
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REACTIONS (INCLUDE IN STRUCTURE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in reaction information into the model structure. Replace colliding 
% parameter names with their prefixed versions (reaction name). Reactions 
% are included into the model structure in the 'reactions' field. This leads 
% to shorter ODE definitions for the specie states, which then contain the 
% reaction names with the correct stoichiometries. The 'reactions' field
% contains three fields: 'name', 'formula', 'notes'
%
% We differentiate the case where:
%   - single compartment model with constant compartment size = 1
%     (nothing is done additionally)
%   - multi compartment model or single compartment model with compartment 
%     size different from one. In this case all species names in the
%     reactions are exchanged agains speciesname/speciesCompartmentSize
%
% Cycle through all reactions
for k1 = 1:length(SBMLmodel.reaction),
    reactionName = SBMLmodel.reaction(k1).id;
    kineticLaw = SBMLmodel.reaction(k1).kineticLaw;
    % check if a kineticLaw is given
    if isempty(kineticLaw),
        errorMsg = sprintf('%s\nNo kinetic law is given for reaction ''%s''\n',errorMsg,reactionName);
    elseif isempty(kineticLaw.formula) || isempty(kineticLaw.math),
        errorMsg = sprintf('%s\nNo kinetic formula law given for reaction ''%s''\n',errorMsg,reactionName);
    end
    % Process the kineticLaw formula by replacing the local parameter names
    % by reactionname_parameternames and replace the mathml expressions
    if ~isempty(kineticLaw),
        if length(kineticLaw.formula) >= length(kineticLaw.math),
            formula = kineticLaw.formula;
        else
            formula = kineticLaw.math;
        end
        if ~isempty(kineticLaw.parameter),
            % SBML 2 can have formula in formula or math. Check which string is
            % longer and use this one (quick and dirty fix - since it has been
            % observer that TranslateSBML can return strange stuff in the math
            % field sometimes)
            for k2 = 1:length(kineticLaw.parameter),
                parameterName = kineticLaw.parameter(k2).id;
                newParameterName = char([double(reactionName) double('_') double(parameterName)]); % fastest strcat
                formula = exchangeStringInString(formula,parameterName,newParameterName);
            end
        end
    else
        formula = '';
    end
    % replace MathML expressions in formula against MATLAB expressions
    [formula] = replaceMathMLexpressions(formula);
    % Include reaction information in the structure
    IQMstructure.reactions(k1).name = reactionName;
    IQMstructure.reactions(k1).formula = formula;
    IQMstructure.reactions(k1).notes = convert2IQMNotes(SBMLmodel.reaction(k1).notes,1);
    IQMstructure.reactions(k1).reversible = SBMLmodel.reaction(k1).reversible;    
    if (SBMLmodel.reaction(k1).isSetFast == 1) && (SBMLmodel.reaction(k1).fast == 1),
        IQMstructure.reactions(k1).fast = 1;    
    else
        IQMstructure.reactions(k1).fast = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS REACTION ODES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add the reaction variable names with correct stoichiometries to the 'ODE' field
% for the correct species. The 'correct' species are only non-boundary and 
% non-constant species for which no rules exist. These are listed in
% the 'speciesODElist'
%
for k1 = 1:length(speciesODElist),
    % get name of species and the index of the corresponding state
    speciesName = speciesODElist(k1).name;
    stateIndex = speciesODElist(k1).stateIndex;
    % cycle through all the reactions and look in which reactions the
    % species is changed
    for k2 = 1:length(IQMstructure.reactions),
        % get reaction name for inclusion in ODEstring
        reactionName = IQMstructure.reactions(k2).name;
        % cycle through the reactants of the current reaction
        % (ordering of reactions in IQMstructure and SBMLmodel are equal)
        for k3 = 1:length(SBMLmodel.reaction(k2).reactant),
            % Handle missing stoichiomety definition in L2V>1 - then set to 1
            if isnan(SBMLmodel.reaction(k2).reactant(k3).stoichiometry),
                SBMLmodel.reaction(k2).reactant(k3).stoichiometry = 1;
            end
            % Handle missing denominator definition in L2V>1 - then set to 1
            if ~isfield(SBMLmodel.reaction(k2).reactant(k3),'denominator'),
                SBMLmodel.reaction(k2).reactant(k3).denominator = 1;
            else
                if isempty(SBMLmodel.reaction(k2).reactant(k3).denominator),
                    SBMLmodel.reaction(k2).reactant(k3).denominator = 1;
                end
            end
            
            reactantName = SBMLmodel.reaction(k2).reactant(k3).species;
            if ~isempty(SBMLmodel.reaction(k2).reactant(k3).stoichiometryMath),
                reactantStoichiometry = SBMLmodel.reaction(k2).reactant(k3).stoichiometryMath; 
                if isstruct(reactantStoichiometry),
                    reactantStoichiometry = reactantStoichiometry.math;
                end                
                isunitystoich = 0;
            else
                reactantStoichiometry_num = abs(double(SBMLmodel.reaction(k2).reactant(k3).stoichiometry) / double(SBMLmodel.reaction(k2).reactant(k3).denominator));
                isunitystoich = reactantStoichiometry_num == 1;
                reactantStoichiometry = num2str(reactantStoichiometry_num);
            end
            % If species name and reactant name are equal then add reaction
            % to species ODE
            if strcmp(reactantName,speciesName),
                % construct the string to add to the ODE of the current reactant
                if isunitystoich == 0,
                    ODEstringAdd = char([double('-') double(reactantStoichiometry) double('*') double(reactionName)]); % fastest strcat
                else
                    ODEstringAdd = char([double('-') double(reactionName)]); % fastest strcat
                end
                % add the ODEstringAdd to the 'ODE' field of this species
                IQMstructure.states(stateIndex).ODE = char([double(IQMstructure.states(stateIndex).ODE) double(ODEstringAdd)]); % fastest strcat
                break;
            end
        end
        % cycle through the products of the current reaction
        % (ordering of reactions in IQMstructure and SBMLmodel are equal)
        for k3 = 1:length(SBMLmodel.reaction(k2).product),
            % Handle missing stoichiomety definition in L2V>1 - then set to 1
            if isnan(SBMLmodel.reaction(k2).product(k3).stoichiometry),
                SBMLmodel.reaction(k2).product(k3).stoichiometry = 1;
            end
            % Handle missing denominator definition in L2V>1 - then set to 1
            if ~isfield(SBMLmodel.reaction(k2).product(k3),'denominator'),
                SBMLmodel.reaction(k2).product(k3).denominator = 1;
            else
                if isempty(SBMLmodel.reaction(k2).product(k3).denominator),
                    SBMLmodel.reaction(k2).product(k3).denominator = 1;
                end
            end
            
            productName = SBMLmodel.reaction(k2).product(k3).species;
            if ~isempty(SBMLmodel.reaction(k2).product(k3).stoichiometryMath),
                productStoichiometry = SBMLmodel.reaction(k2).product(k3).stoichiometryMath; 
                if isstruct(productStoichiometry),
                    productStoichiometry = productStoichiometry.math;
                end
                isunitystoich = 0;
            else
                productStoichiometry_num = abs(double(SBMLmodel.reaction(k2).product(k3).stoichiometry) / double(SBMLmodel.reaction(k2).product(k3).denominator));
                isunitystoich = productStoichiometry_num == 1;
                productStoichiometry = num2str(productStoichiometry_num);
            end
            % If species name and product name are equal then add reaction
            % to species ODE
            if strcmp(productName,speciesName),
                % construct the string to add to the ODE of the current product
                if isunitystoich == 0,
                    ODEstringAdd = char([double('+') double(productStoichiometry) double('*') double(reactionName)]); % fastest strcat
                else
                    ODEstringAdd = char([double('+') double(reactionName)]); % fastest strcat
                end
                % add the ODEstringAdd to the 'ODE' field of this species
                IQMstructure.states(stateIndex).ODE = char([double(IQMstructure.states(stateIndex).ODE) double(ODEstringAdd)]); % fastest strcat
                break;
            end
        end
    end
end
% Now all ODEs should be defined - check if everything is correct!
% Finally cycle through all the states and check if all ODEs have been defined
for k = 1:length(IQMstructure.states),
    stateName = IQMstructure.states(k).name;
    stateODE = IQMstructure.states(k).ODE;
    if isempty(stateODE),
%         % Issue a warning only and set ODE to 0
%         disp(sprintf('No ODE defined for state ''%s''',stateName));
        IQMstructure.states(k).ODE = '0';
    end
end
% After initial construction of the ODEs their units need to be changed (or
% not) since the rate laws are assumed to have amount/time units.
%
% - species defined as 'amount' => no change necessary!
% - species defined as 'concentration' => divide the ode expression by the
%   compartment size the species is in.
%
% We differentiate the case where:
%   - single compartment model with constant compartment size = 1
%     (nothing is done additionally)
%   - multi compartment model or single compartment model with compartment 
%     size different from one. In this case the conversion is done for the
%     concentration species.
%
compartmentddtdonelist = {};
for k1 = 1:length(speciesODElist),
    % get name of species and the index of the corresponding state
    speciesName = speciesODElist(k1).name;
    speciesUnits = speciesODElist(k1).units;
    stateIndex = speciesODElist(k1).stateIndex;
    % only do the conversion if species in concentration units!
    if strcmp(speciesUnits,'concentration'),
        % check if the model has only one compartment with volume 1 and
        % constant volume => if yes then leave ODE as it is.
        convertODE = 1;
%         if length(SBMLmodel.compartment) == 1,
%             % only one compartment
%             if SBMLmodel.compartment(1).size == 1,
%                 % size of compartment is 1
%                 for k = 1:length(IQMstructure.parameters),
%                     % check if compartment size is constant (then it is defined as
%                     % a parameter)
%                     if strcmp(SBMLmodel.compartment(1).id,IQMstructure.parameters(k).name),
%                         % compartment defined as a parameter and thus no species
%                         % names need to be exchanged
%                         convertODE = 0;
%                         % return directly to the main function
%                     end
%                 end
%             end
%         end
        if convertODE,
            % convert the ODE expression
            % get compartment name for species
            compartmentName = getCompartmentNameSpecies(speciesName,SBMLmodel);
            ODEstring = IQMstructure.states(stateIndex).ODE;
            if ~strcmp(ODEstring,'0')
                % construct new ODE string
                newODEstring = char([double('(') double(ODEstring) double(')/') double(compartmentName)]); % fastest strcat 
                % include the new ODE string in structure
                IQMstructure.states(stateIndex).ODE = newODEstring;
            end
            
            
% HANDLE TIME VARYING COMPARTMENT SIZES
            [changingFlag,ddtCompartment] = getCompartmentInfoChanging(compartmentName,IQMstructure,compartmentddtdonelist);
            if changingFlag,
                % append the second derivative term
                IQMstructure.states(stateIndex).ODE = sprintf('%s - %s*ddt_%s/%s',IQMstructure.states(stateIndex).ODE,speciesName,compartmentName,compartmentName);
                
                % add the time derivative of the compartment to variables
                % (if not done already)
                if isempty(strmatchIQM(compartmentName,compartmentddtdonelist,'exact')),
                    IQMstructure.variables(end+1).name = ['ddt_' compartmentName];
                    IQMstructure.variables(end).formula = ddtCompartment;
                    IQMstructure.variables(end).type = 'isParameter';
                    IQMstructure.variables(end).compartment = '';
                    IQMstructure.variables(end).unittype = '';
                    IQMstructure.variables(end).notes = 'Time derivative of compartment';
                    compartmentddtdonelist{end+1} = compartmentName;
                end
            end

            
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EVENTS (INCLUDE IN STRUCTURE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k1 = 1:length(SBMLmodel.event),
    % get the name of the event
    if ~isempty(SBMLmodel.event(k1).id),
        eventName = SBMLmodel.event(k1).id;
    else 
        eventName = sprintf('Event_%d',k1);
    end
    % check event trigger is present
    if isempty(SBMLmodel.event(k1).trigger),
        errorMsg = sprintf('%s\nNo trigger condition defined for event ''%s''.',errorMsg,eventName);
    end
    % get delay and check if it is not given or given
    delay = SBMLmodel.event(k1).delay;
    if isstruct(delay) && ~isempty(delay),
        delay = delay.math;
    end
    if isempty(delay),
        delay = '0';
    end
    delay = replaceMathMLexpressions(delay);
    % fill structure with basic event information
    IQMstructure.events(k1).name = eventName;
    % check the trigger expression and handle delay if present
    trigger = SBMLmodel.event(k1).trigger;
    if isstruct(trigger),
        % Happens in L2 V3 or definitly 4
        trigger = trigger.math;
    end
        

    trigger = removeWhiteSpace(trigger);
    trigger = replaceMathMLexpressions(trigger);
    if ~strcmp(delay,'0'),
        trigger = ['delayIQM(' trigger ',' delay ')'];
    end
    IQMstructure.events(k1).trigger = trigger;
    IQMstructure.events(k1).notes = convert2IQMNotes(SBMLmodel.event(k1).notes,1);
    % fill in data for event assignments
    for k2 = 1:length(SBMLmodel.event(k1).eventAssignment),
        variableName = SBMLmodel.event(k1).eventAssignment(k2).variable;
        variableValue = SBMLmodel.event(k1).eventAssignment(k2).math;
        % check that the name is the name of a state in the model
        nameFound = 0;
        for k3 = 1:length(IQMstructure.states),
            if strcmp(variableName,IQMstructure.states(k3).name),
                nameFound = 1;
                break;
            end
        end
        if nameFound == 0,
            % events on parameters are now allowed
            %errorMsg = sprintf('%s\nThe variable affected by event ''%s'' is not a state.',errorMsg,eventName);
        end
        % check that the variable value is set
        if isempty(variableValue),
            errorMsg = sprintf('%s\nNo math element given for event ''%s''.',errorMsg,eventName);
        end
        % update IQMstructure
        IQMstructure.events(k1).assignment(k2).variable = variableName;
        assignformula = replaceMathMLexpressions(variableValue);
        % if delay not 0 then add delay information
        if ~strcmp(delay,'0'),
            assignformula = ['delayIQM(' assignformula ',' delay ')'];
        end
        IQMstructure.events(k1).assignment(k2).formula = assignformula;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXCHANGE NaN VALUES AGAINST 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:length(IQMstructure.states),
    if isnan(IQMstructure.states(k).initialCondition),
        IQMstructure.states(k).initialCondition = 0;
    end
end
for k=1:length(IQMstructure.parameters),
    if isnan(IQMstructure.parameters(k).value),
        IQMstructure.parameters(k).value = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORT ALGEBRAIC RULES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare the import
% 1)need to find all species, compartments, and (global) parameters that are
% NOT constant. These are the candidates for the determination by an
% algebraic rule (in this order of priority)
nonConstantSpecies = {};
nonConstantParameters = {};
nonConstantCompartments = {};
for k=1:length(SBMLmodel.species),
    if SBMLmodel.species(k).constant == 0,
        nonConstantSpecies{end+1} = SBMLmodel.species(k).id;
    end
end
for k=1:length(SBMLmodel.parameter),
    if SBMLmodel.parameter(k).constant == 0,
        nonConstantParameters{end+1} = SBMLmodel.parameter(k).id;
    end
end
for k=1:length(SBMLmodel.compartment),
    if SBMLmodel.compartment(k).constant == 0,
        nonConstantCompartments{end+1} = SBMLmodel.compartment(k).id;
    end
end
% 2) now we need to make a list of all states determined by ODEs and
% variables determined by formulas. 
allStates = {};
for k=1:length(IQMstructure.states),
    if IQMstructure.states(k).ODE ~= '0',
        % Only non constant states (would also require the check of event 
        % assignments ... but in my opinion, SBML should simply define the
        % variable which is affected by an algebraic rule ... so much
        % cleaner!
        allStates{end+1} = IQMstructure.states(k).name;
    end
end
allVariables = {};
for k=1:length(IQMstructure.variables),
    if ~isnan(str2double(IQMstructure.variables(k).formula)),
        % If variable has not a numeric assignment ... 
        allVariables{end+1} = IQMstructure.variables(k).name;
    end
end
% 3) The names of the states and variables need to be deleted from the
% nonConstantSpecies, nonConstantCompartments, and nonConstantParameters,
% because only the remaining ones can be determined by algebraic rules
possibleSpecies = setdiff(nonConstantSpecies,allStates);
possibleSpecies = setdiff(possibleSpecies,allVariables);
possibleParameters = setdiff(nonConstantParameters,allStates);
possibleParameters = setdiff(possibleParameters,allVariables);
possibleCompartments = setdiff(nonConstantCompartments,allStates);
possibleCompartments = setdiff(possibleCompartments,allVariables);
% 4) now we check if at least as many possibilities as there are ARs
% (otherwise certainly overdetermined).
ARindices = strmatchIQM('SBML_ALGEBRAIC_RULE',{SBMLmodel.rule.typecode},'exact');
nrARs = length(ARindices);
if length(possibleSpecies) + length(possibleParameters) + length(possibleCompartments) < nrARs,
    error('The model is certainly overdetermined. To many algebraic rules!');
end
% 5) Now check which of the possible elements appear in the ARs. If at
% least a number of nrARs do appear there its fine ... if not put out a
% warning and select whatever.
ARs = SBMLmodel.rule(ARindices);
ReallyPossibleSpecies = {};
ReallyPossibleParameters = {};
ReallyPossibleCompartments = {};
for k=1:nrARs,
    species = regexp(ARs(k).formula,possibleSpecies,'match');
    for k2=1:length(species),
        if ~isempty(species{k2}),
            ReallyPossibleSpecies{end+1} = species{k2}{1};
        end
    end
    parameters = regexp(ARs(k).formula,possibleParameters,'match');
    for k2=1:length(parameters),
        if ~isempty(parameters{k2}),
            ReallyPossibleParameters{end+1} = parameters{k2}{1};
        end
    end
    compartments = regexp(ARs(k).formula,possibleCompartments,'match');
    for k2=1:length(compartments),
        if ~isempty(compartments{k2}),
            ReallyPossibleCompartments{end+1} = compartments{k2}{1};
        end
    end
end
% make unique (it might not be)
ReallyPossibleSpecies = unique(ReallyPossibleSpecies);
ReallyPossibleParameters = unique(ReallyPossibleParameters);
ReallyPossibleCompartments = unique(ReallyPossibleCompartments);
nrPossible = length(ReallyPossibleSpecies) + length(ReallyPossibleParameters) + length(ReallyPossibleCompartments);
if nrPossible < nrARs,
    error('The model is overdetermined or at least it is unclear which variables are to be determined using the algebraic rules.');
end
% get the indices in the SBML model for the really possible ones
indexspecies = [];
indexparameters = [];
indexcompartments = [];
for k=1:length(ReallyPossibleSpecies),
    indexspecies(end+1) = strmatchIQM(ReallyPossibleSpecies{k},{SBMLmodel.species.id},'exact');
end
for k=1:length(ReallyPossibleParameters),
    indexparameters(end+1) = strmatchIQM(ReallyPossibleParameters{k},{SBMLmodel.parameter.id},'exact');
end
for k=1:length(ReallyPossibleCompartments),
    indexcompartments(end+1) = strmatchIQM(ReallyPossibleCompartments{k},{SBMLmodel.compartment.id},'exact');
end
% 6) Add the algebraic rules and delete the corresponding parameters from
% the IQMmodel structure
for k=1:length(ARindices),
    rule = SBMLmodel.rule(ARindices);
    % determine the name (it needs to be one of the elements in
    % indexspecies, indexparameters, or indexcompartments
    if ~isempty(indexspecies),      
        index = indexspecies(1);
        name = SBMLmodel.species(index).id;
        type = 'isSpecie';
        compartment = SBMLmodel.species(index).compartment;
        if SBMLmodel.species(index).isSetInitialAmount,
            if SBMLmodel.species(index).hasOnlySubstanceUnits,
                unittype = 'amount';        
            else
                unittype = 'concentration';
                amount2concentration(end+1).species = name;
                amount2concentration(end).compartment = compartment;                
            end
        else
            unittype = 'concentration';
        end
        % remove the used species
        indexspecies(1) = [];             
    elseif ~isempty(indexparameters),
        index = indexparameters(1);
        name = SBMLmodel.parameter(index).id;
        type = 'isParameter';
        compartment = '';
        unittype = '';
        % remove the used parameter
        indexparameters(1) = [];      
    elseif ~isempty(indexcompartments),
        index = indexcompartments(1);
        name = SBMLmodel.compartment(index).id;
        type = 'isCompartment';
        compartment = SBMLmodel.compartment(index).outside;
        unittype = '';        
        % remove the used compartment
        indexcompartments(1) = [];      
    else
        error('Should not get to this :)');
    end
    % get initial condition (can be state, parameter, or variable .... YEEEESSS!)
    % and delete the corresponding component from the model (now determined
    % by algebraic rule).
    allStates = {IQMstructure.states.name};    
    allParams = {IQMstructure.parameters.name};
    allVars = {IQMstructure.variables.name};
    indexstate = strmatchIQM(name,allStates,'exact');
    indexparam = strmatchIQM(name,allParams,'exact');
    indexvar = strmatchIQM(name,allVars,'exact');
    if ~isempty(indexstate),
        initialCondition = IQMstructure.states(indexstate).initialCondition;    
        IQMstructure.states(indexstate) = [];
    elseif ~isempty(indexparam),
        initialCondition = IQMstructure.parameters(indexparam).value;
        IQMstructure.parameters(indexparam) = [];
    elseif ~isempty(indexvar),
        initialCondition = str2double(IQMstructure.variables(indexvar).formula);    
        IQMstructure.variables(indexvar) = [];
    else
        error('Should not get to this :) ... 2');
    end
    % add algebraic rule
    IQMstructure.algebraic(k).name = name;
    IQMstructure.algebraic(k).formula = rule(k).formula;
    IQMstructure.algebraic(k).initialCondition = initialCondition;
    IQMstructure.algebraic(k).type = type;
    IQMstructure.algebraic(k).compartment = compartment;
    IQMstructure.algebraic(k).unittype = unittype;
    IQMstructure.algebraic(k).notes = rule(k).notes;
end
% 7) Its DONE!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK ALL COMPONENTNAMES AND REMOVE TRAILING "_"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get names to exchange
oldnames = {};
newnames = {};
% states
for k=1:length(IQMstructure.states),
    [oldnames,newnames,changed,IQMstructure.states(k).name] = removeUnderscores(oldnames,newnames,IQMstructure.states(k).name);    
end
% algebraic
for k=1:length(IQMstructure.algebraic),
    [oldnames,newnames,changed,IQMstructure.algebraic(k).name] = removeUnderscores(oldnames,newnames,IQMstructure.algebraic(k).name);    
end
% parameters
for k=1:length(IQMstructure.parameters),
    [oldnames,newnames,changed,IQMstructure.parameters(k).name] = removeUnderscores(oldnames,newnames,IQMstructure.parameters(k).name);
end
% variables
for k=1:length(IQMstructure.variables),
    [oldnames,newnames,changed,IQMstructure.variables(k).name] = removeUnderscores(oldnames,newnames,IQMstructure.variables(k).name);
end
% reactions
for k=1:length(IQMstructure.reactions),
    [oldnames,newnames,changed,IQMstructure.reactions(k).name] = removeUnderscores(oldnames,newnames,IQMstructure.reactions(k).name);   
end
% functions
for k=1:length(IQMstructure.functions),
    [oldnames,newnames,changed,IQMstructure.functions(k).name] = removeUnderscores(oldnames,newnames,IQMstructure.functions(k).name);   
end
% ok, if names changed then handle formulas
if ~isempty(oldnames),
    % states
    for k=1:length(IQMstructure.states),
        IQMstructure.states(k).ODE = regexprep(IQMstructure.states(k).ODE,oldnames,newnames);
        IQMstructure.states(k).compartment = regexprep(IQMstructure.states(k).compartment,oldnames,newnames);
    end
    % algebraic
    for k=1:length(IQMstructure.algebraic),
        IQMstructure.algebraic(k).formula = regexprep(IQMstructure.algebraic(k).formula,oldnames,newnames);    
        IQMstructure.algebraic(k).compartment = regexprep(IQMstructure.algebraic(k).compartment,oldnames,newnames);
    end
    % parameters
    for k=1:length(IQMstructure.parameters),
        IQMstructure.parameters(k).compartment = regexprep(IQMstructure.parameters(k).compartment,oldnames,newnames);
    end    
    % variables
    for k=1:length(IQMstructure.variables),
        IQMstructure.variables(k).formula = regexprep(IQMstructure.variables(k).formula,oldnames,newnames);
        IQMstructure.variables(k).compartment = regexprep(IQMstructure.variables(k).compartment,oldnames,newnames);
    end
    % reactions
    for k=1:length(IQMstructure.reactions),
        IQMstructure.reactions(k).formula = regexprep(IQMstructure.reactions(k).formula,oldnames,newnames);
    end
    % functions
    for k=1:length(IQMstructure.functions),
        IQMstructure.functions(k).formula = regexprep(IQMstructure.functions(k).formula,oldnames,newnames);
    end
    % events
    for k=1:length(IQMstructure.events),
        IQMstructure.events(k).trigger = regexprep(IQMstructure.events(k).trigger,oldnames,newnames);
        for k2=1:length(IQMstructure.events(k).assignment),
            IQMstructure.events(k).assignment(k2).variable = regexprep(IQMstructure.events(k).assignment(k2).variable,oldnames,newnames);
            IQMstructure.events(k).assignment(k2).formula = regexprep(IQMstructure.events(k).assignment(k2).formula,oldnames,newnames);
        end
    end
end

%% Here handle initial assignments
% Idea: 
% - get the name of the element that has an initial assignment.
% - for states implement as initial condition
% - for variables as new formula (valid until time=0?)
% - for parameters: make variables out of parameters and handle as variables

if ~isfield(SBMLmodel,'initialAssignment'),
    % Add field if not present
    SBMLmodel.initialAssignment = [];
end

% Cycle through all initial assignments
for k=1:length(SBMLmodel.initialAssignment),
    symbol = SBMLmodel.initialAssignment(k).symbol;
    math = replaceMathMLexpressions(SBMLmodel.initialAssignment(k).math);
    % Adjust math if species and in amount 
    ix = strmatchIQM(symbol,speciesAmount_SAVEnames,'exact');
    if ~isempty(ix),
        math = ['(' math ')*' speciesAmount_SAVEcompartments{ix}];
    end
    % Handle states
    ix = strmatchIQM(symbol,{IQMstructure.states.name},'exact');
    if ~isempty(ix),
        IQMstructure.states(ix).initialCondition = math;
    end
    % Handle variables
    ix = strmatchIQM(symbol,{IQMstructure.variables.name},'exact');
    if ~isempty(ix),
        IQMstructure.variables(ix).formula = math;
    end    
    % Handle parameters
    ix = strmatchIQM(symbol,{IQMstructure.parameters.name},'exact');
    if ~isempty(ix),
        paramsave = IQMstructure.parameters(ix);
        % Remove parameter from model
        IQMstructure.parameters(ix) = [];
        % Create variable
        variablenew = [];
        variablenew.name        = paramsave.name;
        variablenew.formula     = math;
        variablenew.type        = paramsave.type;
        variablenew.compartment = paramsave.compartment;
        variablenew.unittype    = paramsave.unittype;
        variablenew.notes       = paramsave.notes;
        % Add variable
        IQMstructure.variables(end+1) = variablenew;
    end    
    % Handle reactions
    ix = strmatchIQM(symbol,{IQMstructure.reactions.name},'exact');
    if ~isempty(ix),
        error('Initial Assignments on Reactions do not make sense.');
    end    
end

%% Handle species that are in amount but in expressions in concentration
for k0=1:length(speciesAmount_SAVEnames),
    name = speciesAmount_SAVEnames{k0};
    comp = speciesAmount_SAVEcompartments{k0};
    replacename = ['(' name '/' comp ')'];
    % In ODEs
    for k=1:length(IQMstructure.states),
        IQMstructure.states(k).ODE = regexprep(IQMstructure.states(k).ODE,['\<' name '\>'],replacename);
    end
    % In initial conditions of ODEs
    for k=1:length(IQMstructure.states),
        if ~isnumeric(IQMstructure.states(k).initialCondition),
            IQMstructure.states(k).initialCondition = regexprep(IQMstructure.states(k).initialCondition,['\<' name '\>'],replacename);
        end
    end
    % In variables
    for k=1:length(IQMstructure.variables),
        IQMstructure.variables(k).formula = regexprep(IQMstructure.variables(k).formula,['\<' name '\>'],replacename);
    end
    % In reactions
    for k=1:length(IQMstructure.reactions),
        IQMstructure.reactions(k).formula = regexprep(IQMstructure.reactions(k).formula,['\<' name '\>'],replacename);
    end
    % In event triggers
    for k=1:length(IQMstructure.events),
        IQMstructure.events(k).trigger = regexprep(IQMstructure.events(k).trigger,['\<' name '\>'],replacename);
        % Adjust event assignments
        for k2=1:length(IQMstructure.events(k).assignment),
            if strcmp(IQMstructure.events(k).assignment(k2).variable,name),
                IQMstructure.events(k).assignment(k2).formula = ['(' IQMstructure.events(k).assignment(k2).formula ')*' comp];
            end
        end
    end    
end
    
%% Handle parameters that are changed by an event ... call these parameters differently and create a variable
% with the parameter name ... important to get the test suite to work for few models
parChangeNames = {};
for k=1:length(IQMstructure.events),
    for k2=1:length(IQMstructure.events(k).assignment),
        name = IQMstructure.events(k).assignment(k2).variable;
        % check if parameter
        ix = strmatchIQM(name,{IQMstructure.parameters.name},'exact');
        if ~isempty(ix),
            parChangeNames{end+1} = name;
        end
    end
end
parChangeNames = unique(parChangeNames);

for k=1:length(parChangeNames),
    name                    = parChangeNames{k};
    % Get param info and change its name to "name"_par
    ix                      = strmatchIQM(name,{IQMstructure.parameters.name},'exact');
    parstruct               = IQMstructure.parameters(ix);
    IQMstructure.parameters(ix).name = [name '_par'];
    % Create a variabe "name" = "name"_par
    variablenew             = [];
    variablenew.name        = name;
    variablenew.formula     = [name '_par'];
    variablenew.type        = parstruct.type;
    variablenew.compartment = parstruct.compartment;
    variablenew.unittype    = parstruct.unittype;
    variablenew.notes       = parstruct.notes;
    % Add variable as first into the model structure
    allVariables            = IQMstructure.variables;    
    if isempty(allVariables),
        variableall         = variablenew;
    else
        % Create a dummy first and shift the present ones to start at index 2
        variableall         = allVariables([1 1:length(allVariables)]);
        % Add the new variable
        variableall(1)      = variablenew;
    end
    IQMstructure.variables  = variableall;
    % Change event assignment par names
    for k2=1:length(IQMstructure.events),
        for k3=1:length(IQMstructure.events(k2).assignment),
            IQMstructure.events(k2).assignment(k3).variable = regexprep(IQMstructure.events(k2).assignment(k3).variable, ['\<' name '\>'],[name '_par']);
        end
    end
end

%% return from main function
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLE TIME VARYING COMPARTMENTS ... get changing flag and the
% time-derivative of the compartment size.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [changingFlag,ddtCompartment] = getCompartmentInfoChanging(compartmentName,IQMstructure,compartmentddtdonelist)
changingFlag = 0; 
ddtCompartment = '';
% go through states and variables to check if the compartment name appears
% there ... then just do it
indexState = strmatchIQM(compartmentName,{IQMstructure.states.name},'exact');
indexVariable = strmatchIQM(compartmentName,{IQMstructure.variables.name},'exact');
indexAlgebraic = strmatchIQM(compartmentName,{IQMstructure.algebraic.name},'exact');
if ~isempty(indexState),
    changingFlag = 1; 
    ddtCompartment = IQMstructure.states(indexState).ODE;
elseif ~isempty(indexVariable),
    changingFlag = 1;
    Compartment = IQMstructure.variables(indexVariable).formula;
    % ddtCompartment needs only to be determined once for each compartment.
    % So not per species in this compartment. We make it simple here (since
    % it mainly matters when determining and adding it manually):
    if isempty(strmatchIQM(compartmentName,compartmentddtdonelist,'exact')),
        if isSymbolicpresentIQM(),
            ddtCompartment = char(diff(Compartment,'time'));
        else
            disp('For the handling of time varying compartments the time derivative of the');
            disp('compartment size needs to be determined. Since the symbolic toolbox is not');
            disp('present you need to do that manually.');
            disp(' ');
            disp(sprintf('C(time) = %s',Compartment));
            disp(' ');
            ddtCompartment = input('d/dtime(C(time)) = ','s');
        end
    end
elseif ~isempty(indexAlgebraic),
    error(sprintf('Size of compartment ''%s'' defined by an algebraic rule! Cant determine\nthe analytic rate of change of the size.',compartmentName));
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELP FUNCTION TO REMOVE __ at the start of element names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [oldnames,newnames,changed,newname] = removeUnderscores(oldnames,newnames,name)
changed = 0;
newname = regexprep(name,'^([_]?)','underscore_');
if length(name) ~= length(newname),
    changed = 1;
    oldnames{end+1} = sprintf('\\<%s\\>',name);
    newnames{end+1} = newname;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND RULES FOR A CERTAIN VARIABLE NAME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cycle through all the rules and check if the given variable name
% (species or parameter or compartment) appears as left hand side (LHS) in a rule
%
% numberRulesVariable: number of rules the name appears in as LHS
% lastRuleIndex: the index of the rule the name was last detected
% lastRuleType: 'scalar' or 'rate'
% lastRuleFormula: the RHS expression for the variable
%
% We assume that the name of the variable can appear in the
% 'variable' field, in the 'species' field, in the 'compartment' field,
% or in the 'name' field of the 'rule' field. 
% Some SBML example models had species in the 'variable' field!
%
function [numberRulesVariable,lastRuleIndex,lastRuleType,lastRuleFormula,errorMsg] = getRules(variableName,SBMLmodel,errorMsg)
% Initialize return values
numberRulesVariable = 0;
lastRuleIndex = [];
lastRuleType = '';
lastRuleFormula = '';
for k = 1:length(SBMLmodel.rule),
    % DO NOT PROCESS ALGEBRAIC RULES HERE ... WE DO THAT AT A DIFFERENT PLACE
    if isempty(strfind(SBMLmodel.rule(k).typecode,'ALGEBRAIC')),
        if strcmp(variableName,SBMLmodel.rule(k).variable) || strcmp(variableName,SBMLmodel.rule(k).species) || strcmp(variableName,SBMLmodel.rule(k).name) || strcmp(variableName,SBMLmodel.rule(k).compartment),
            % The name was found as LHS in a rule
            numberRulesVariable = numberRulesVariable+1;
            lastRuleIndex = k;
            if isempty(strfind(SBMLmodel.rule(k).typecode,'RATE')),
                % Rule is a scalar rule
                lastRuleType = 'scalar';
            else
                % Rule is a rate rule
                lastRuleType = 'rate';
            end
            % get the RHS expression of the rule
            lastRuleFormula = SBMLmodel.rule(k).formula;
        end
    end
end
% replace MathML expressions against MATLAB expressions
[lastRuleFormula] = replaceMathMLexpressions(lastRuleFormula);
% return
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK IF A CERTAIN SPECIES APPEARS AS PRODUCT AND/OR REACTANT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In the case that a non-boundary species is defined by a rule
% we need to check if it is also present as product or reactant in any
% reaction. In this case the SBML model is not correctly defined and
% reactionSpeciesPresent = 1
function [reactionSpeciesPresent] = checkReactionsForSpecies(speciesName,SBMLmodel)
% Initialize return value
reactionSpeciesPresent = 0;
reactionComponents = {};
for k = 1:length(SBMLmodel.reaction),
    reactionComponents = {reactionComponents{:},SBMLmodel.reaction(k).reactant.species,SBMLmodel.reaction(k).product.species};
end
if ~isempty(strmatchIQM(speciesName,reactionComponents,'exact')),
    reactionSpeciesPresent = 1;
end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DELETE WHITESPACES IN STRINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Useful for taking away whitespaces in kineticLaw formulas, as
% seen in some example models
function [outputString] = removeWhiteSpace(inputString)
outputString = strrep(inputString,' ','');
% return
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET COMPARTMENT NAME FOR SPECIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [compartmentName] = getCompartmentNameSpecies(species,SBMLmodel)
compartmentName = SBMLmodel.species(strmatchIQM(species,{SBMLmodel.species.id},'exact')).compartment;
% return
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXCHANGE STRING IN STRING - WITH CHECK OF ALLOWED CHARACTERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find occurrences of stringOld in fullString and exchange it with
% stringNew. 
function [processedString] = exchangeStringInString(fullString,stringOld,stringNew);
% really simple!!!
exprString = char([double('\<') double(stringOld) double('\>')]);
processedString = regexprep(fullString, exprString, stringNew);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXCHANGE MATHML FUNCTION NAMES IN ALL FORMULAS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MathML functions can have different names than MATLAB function names
% we need to exchange them in the formula strings. The term 'MathML
% expression' relates to what is returned by TranslateSBML!
function [newFormula] = replaceMathMLexpressions(formula)
% MathML expressions that need simple exchange with corresponding MATLAB
% expressions:
MathMLexpressions = {'\<multiply\>','\<delay\>','\<piecewise\>','\<and\>','\<or\>','\<arccos\>','\<arcsin\>','\<arctan\>','\<ceiling\>','\<ln\>',...
    '\<pow\>','\<arccos\>','\<arccosh\>','\<arccot\>','\<arccoth\>','\<arccsc\>','\<arccsch\>',...
    '\<arcsec\>','\<arcsech\>','\<arcsin\>','\<arcsinh\>','\<arctan\>','\<arctanh\>',...
    '\<exponentiale','\<geq\>','\<leq\>','\<xor\>','\<multiply\>'};
MATLABexpressions = {'multiplyIQM','delayIQM','piecewiseIQM','andIQM','orIQM','acos','asin','atan','ceil','log', ...
    'power','acos','acosh','acot','acoth','acsc','acsch',...
    'asec','asech','asin','asinh','atan','atanh'...
    'exp(1)','ge','le','xorIQM','multiplyIQM'};
% find indices
newFormula = regexprep(formula,MathMLexpressions,MATLABexpressions);
% replace the time_symbol with 'time' if a time_symbol is defined and not
% 'time'.
newFormula = exchangeTimeSymbol(newFormula);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exchange time symbol if needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [expression] = exchangeTimeSymbol(expression)
global SBMLtimesymbol
if ~isempty(SBMLtimesymbol) && ~strcmp(SBMLtimesymbol,'time'),
    expression = regexprep(expression,strcat('\<',SBMLtimesymbol,'\>'),'time');
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% regexprep command doing the replacement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [formula] = exchangepowerexp(formula)
global changeFlag
oldformula = formula;
formula = regexprep(['#' formula],'([\W]+)',' $1 ');
formula = regexprep(formula,'[\s]power[\s]*\(([^,]+),([^,]+)\)','($1)^($2)');
formula = regexprep(formula,'\s','');
formula = formula(2:end);
if ~strcmp(oldformula,formula),
    changeFlag = 1;
end
return