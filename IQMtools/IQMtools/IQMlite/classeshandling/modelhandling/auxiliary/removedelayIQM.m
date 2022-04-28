function [modelnew] = removedelayIQM(model)
% removedelayIQM: removes the delay function from the model.
% This is for example necessary for the determination of the 
% steady-state. Thus this function is called, e.g., from IQMsteadystate
% in the case that delay functions have been detected in the model.
%
% Output Arguments:
% =================
% modelnew: model with removed delay function

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>

if ~isIQMmodel(model),
    error('Given model is not an IQMmodel.');
end

% get model structure
ms = struct(model);

% check all ODEs
for k=1:length(ms.states),
    ms.states(k).ODE = removeDelay(ms.states(k).ODE);
end
    
% check all variables
for k=1:length(ms.variables),
    ms.variables(k).formula = removeDelay(ms.variables(k).formula);
end

% check all reactions
for k=1:length(ms.reactions),
    ms.reactions(k).formula = removeDelay(ms.reactions(k).formula);
end

% check all event triggers
for k=1:length(ms.events),
    ms.events(k).trigger = removeDelay(ms.events(k).trigger);
end

% check all event assignments
for k=1:length(ms.events),
    for k2=1:length(ms.events(k).assignment),
        ms.events(k).assignment(k2).formula = removeDelay(ms.events(k).assignment(k2).formula);
    end
end

% check all functions
for k=1:length(ms.functions),
    ms.functions(k).formula = removeDelay(ms.functions(k).formula);
end

% check MATLAB functions
ms.functionsMATLAB = removeDelay(ms.functionsMATLAB);

% return new model
modelnew = IQMmodel(ms);
return

% remove the delayIQM function
function [formula] = removeDelay(formula)
% first remove the second input argument to delayIQM
count = 1;
while 1,
    index = strfind(formula,'delayIQM(');
    if length(index) < count,
        break;
    end
    indexstart = index(count)+length('delayIQM(');
    indexend = indexstart;
    % search the end of the delay argument definition
    parOpen = 1;
    while parOpen ~= 0,
        if formula(indexend) == '(',
            parOpen = parOpen + 1;
        elseif formula(indexend) == ')',
            parOpen = parOpen - 1;
        end
        indexend = indexend + 1;
    end
    argument = formula(indexstart:indexend-2);
    % we only need the first argument
    terms = explodePCIQM(argument,',');
    argument = terms{1};
    % reconstruct formula
    firstpart = formula(1:indexstart-1);
    lastpart = formula(indexend-1:end);
    middlepart = argument;
    formula = char([double(firstpart) double(middlepart) double(lastpart)]);
    % increase counters
    count = count + 1;
end
% finally remove the delayIQM call
formula = regexprep(formula,'\<delayIQM\>','');
return
