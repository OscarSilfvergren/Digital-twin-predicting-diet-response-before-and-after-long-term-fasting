function [] = createEventFunctionIQM(iqm,filename,varargin)
% createEventFunctionIQM: writes the event handling function that is
% necessary to pass to the intergrator in order to be able to deal with
% discrete state events when doing simulations. This function is called is
% called by the function IQMcreateODEfile in the case that the event flag is
% set and at least one event is present in the model.
%
% USAGE:
% ======
% [] = createEventFunctionIQM(iqm,filename)
%
% iqm: IQMmodel  (ODE file model description is not possible to use)
% filename: the filename to which wo write the model

if nargin == 3,
    delaybase = varargin{1};
else
    delaybase = '';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK IF IQMmodel (not really needed but safer)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isIQMmodel(iqm),
    error('Function only defined for IQMmodels.');
end
iqmstruct = struct(iqm);

[PATHSTR,functionName,EXT] = fileparts(filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EVENT TRIGGER DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
triggerVector = '';
for k = 1:length(iqmstruct.events),
    delayname = [delaybase '_eventtrigger_' sprintf('%d',k)];    
    trigger = processFormulaIQM(iqmstruct.events(k).trigger,delayname);
    triggerString = sprintf('double(%s)-0.5',trigger);
    % add trigger to vector
    triggerVector = sprintf('%s%s\n',triggerVector,triggerString);
end
triggerVector = triggerVector(1:end-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPEN FILE FOR WRITING AND WRITE HEADER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(filename,'w');
fprintf(fid,'function [value,isterminal,direction] = %s(time_local,statevector,varargin)\n',functionName);
fprintf(fid,'%% This function is passed to the MATLAB integrator in order to detect\n%% events that are present in the model.\n\n');
fprintf(fid,'global time\n');
fprintf(fid,'time = time_local;\n\n');

fprintf(fid,'parameterValuesNew = [];\n');
fprintf(fid,'if nargin == 3,\n');
fprintf(fid,'    parameterValuesNew = varargin{1};\n');
fprintf(fid,'end\n\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITE THE MODEL STUFF TO INITIALIZE all eventually needed data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'%% MODEL DATA\n');
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
% PROCESS STATES
stateNames = IQMstates(iqm);
for k = 1:length(stateNames),
   fprintf(fid,'%s = statevector(%d);\n',stateNames{k},k);
end
% ALGEBRAIC VARIABLES
if ~isempty(iqmstruct.algebraic),
    didWriteStart = 0;
    algebraicNames = {iqmstruct.algebraic.name};
    offset = length(stateNames);
    offset2 = 1;
    for k = 1:length(algebraicNames),
        if didWriteStart == 0 && ~isempty(algebraicNames{k}),
            fprintf(fid,'%% ALGEBRAIC VARIABLES\n');
            didWriteStart = 1;
        end
        if ~isempty(algebraicNames{k}),
            fprintf(fid,'%s = statevector(%d);\n',algebraicNames{k},offset+offset2);
            offset2 = offset2+1;
        end
    end
end


% PROCESS PARAMETERS
[parameterNames,parameterValues] = IQMparameters(iqm);
fprintf(fid,'if isempty(parameterValuesNew),\n');
for k = 1:length(parameterNames),
    fprintf(fid,'\t%s = %g;\n',parameterNames{k},parameterValues(k));
end
fprintf(fid,'else\n'); 
for k = 1:length(parameterNames),
    fprintf(fid,'\t%s = parameterValuesNew(%d);\n',parameterNames{k},k);
end
fprintf(fid,'end\n');
% PROCESS VARIABLES
[variableNames,variableFormulas] = IQMvariables(iqm);
for k = 1:length(variableNames),
   delayname = [delaybase '_var_' sprintf('%d',k)];        
   fprintf(fid,'%s = %s;\n',variableNames{k},processFormulaIQM(variableFormulas{k},delayname));
end 
% PROCESS REACTIONS
[reactionNames,reactionFormulas] = IQMreactions(iqm);
for k = 1:length(reactionNames),
   delayname = [delaybase '_reac_' sprintf('%d',k)];    
   fprintf(fid,'%s = %s;\n',reactionNames{k},processFormulaIQM(reactionFormulas{k},delayname));
end
% PROCESS EVENT ASSIGNMENTS
for k = 1:length(iqmstruct.events),
    for k2 = 1:length(iqmstruct.events(k).assignment),
        delayname = [delaybase '_eventassign_' sprintf('%d',k) '_' sprintf('%d',k2)];
        fprintf(fid,'eventassign_%d_%d = %s;\n',k,k2,processFormulaIQM(iqmstruct.events(k).assignment(k2).formula,delayname));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITE THE EVENT DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'%% EVENT DATA\n');
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'value = [%s]'';\n',triggerVector);
fprintf(fid,'isterminal = ones(1,%d);\n',length(iqmstruct.events));
fprintf(fid,'direction = ones(1,%d);\n',length(iqmstruct.events));
fprintf(fid,'\n');

% return
fprintf(fid,'return\n');
fprintf(fid,'\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Continue writing model stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITE MODEL FUNCTIONS
[functionNames,functionFormulas,functionArguments] = IQMfunctions(iqm);
for k = 1:length(functionNames),
    fprintf(fid,'function [result] = %s(%s)\n',functionNames{k},functionArguments{k});
    fprintf(fid,'global time\n');
    delayname = [delaybase '_func_' sprintf('%d',k)];    
    fprintf(fid,'result = %s;\n',processFormulaIQM(functionFormulas{k},delayname));
    fprintf(fid,'return\n');
    fprintf(fid,'\n');
end
% WRITE THE MATLAB FUNCTIONS
functionsMATLAB = IQMfunctionsMATLAB(iqm);
if ~isempty(functionsMATLAB),
    delayname = [delaybase '_funcmatlab'];    
    fprintf(fid,'%s',processFormulaIQM(functionsMATLAB,delayname));
    fprintf(fid,'\n');
end
% CLOSE FILE
fclose(fid);
return