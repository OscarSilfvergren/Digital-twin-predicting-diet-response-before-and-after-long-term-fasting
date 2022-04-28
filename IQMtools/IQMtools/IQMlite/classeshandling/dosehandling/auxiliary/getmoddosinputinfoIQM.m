function [inputinfo] = getmoddosinputinfoIQM(model,dosing)
% getmoddosinputinfoIQM: Checks the availability of dosing input
% definitions used in the model in the dosing object and returns a
% structure similar to the "input" field structure of an IQMmodel, augmented
% with the dosing information, defined in the IQMdosing object.
%
% If the dosing object does not contain all inputs, required by the model,
% a warning is displayed and the undefined inputs are deleted from the
% "inputinfo" structure. This allows inputs in an IQMmodel to be kept
% undefined.
%
% USAGE:
% ======
% [inputinfo] = getmoddosinputinfoIQM(model,dosing) 
%
% model: IQMmodel
% dose: IQMdosing object
%
% Output Arguments:
% =================
% inputinfo: Same structure as defined in the IQMmodel "input" field.
%   Augmented by the dosing input information, stroed in the dosing object.

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isIQMmodel(model),
    error('First input argument is not an IQMmodel.');
end
if ~isIQMdosing(dosing),
    error('Second input argument is not an IQMdosing object.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get model and dosing structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ms = struct(model);
ds = struct(dosing);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if model contains inputs.
% If not => error!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
if isempty(ms.inputs),
    error('The model does not contain any inputs. Merging with dosing objects does not make sense.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the inputinfo structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputinfo = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run through all inputs, defined in the model 
% and process the information in the dosing
% object accordingly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minputs = ms.inputs;
dinputs = ds.inputs;
for k=1:length(minputs),
    % First get the info stored in the model
    minfo = minputs(k);
    % Check if input exists in the dosing object
    index = strmatchIQM(minfo.name,{dinputs.name},'exact');
    if ~isempty(index),
        % Input does exist in the dosing object => Now add the complete
        % input information to the inputinfo structure 
        % Get the info for this input stored in the dosing object
        dinfo = dinputs(index);
        % Add dosing schedule information to the minfo structure
        minfo.type = dinfo.type;
        minfo.time = dinfo.time;
        minfo.Tlag = dinfo.Tlag;
        minfo.D = dinfo.D;
        minfo.parameters = dinfo.parameters;
        minfo.notes = dinfo.notes;
        % Check if Tlag to be estimated
        minfo.TlagNotes = dinfo.TlagNotes;
        % Add minfo structure to the inputinfo output structure
        if isempty(inputinfo),
            inputinfo = minfo;
        else
            inputinfo(end+1) = minfo;
        end
    else
        % Input does not exist in the dosing object => warning to user 
        % and do not put it into the inputinfo structure
%         disp(sprintf('The input ''%s'', defined in the model, is not defined in the IQMdosing object.',minfo.name));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if any inputs of the models have been 
% found in the dosing object. If not => error!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
if isempty(inputinfo),
    error('No inputs defined in the model could be found in the dosing object. This should be wrong ...')        
end
    