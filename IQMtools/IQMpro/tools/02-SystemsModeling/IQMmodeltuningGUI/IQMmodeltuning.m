function [modeltuned] = IQMmodeltuning(model,data,varargin)
% IQMmodeltuning: Allows to compare and tune a model to one or more sets
% of data. Here: No experiment description is required, the model is
% assumed to already include the experimental settings.
% If no data is available, a time vector can be provided. The model can
% then be tuned. For comparison the simulation result of the nominal model
% is also shown.
%
% USAGE:
% ======
% modeltuned = IQMmodeltuning(model,data)
% modeltuned = IQMmodeltuning(model,data,options)
% modeltuned = IQMmodeltuning(model,time)
%
% model: IQMmodel to tune
% data:  Single IQMmeasurement object or cell-array with IQMmeasurement
%        objects to which to fit the model
% time:  timevector to use for model simulation
% options: still unused
%
% DEFAULT VALUES:
% ===============
%
% Output Arguments:
% =================
% modeltuned: The tuned model with changed parameters

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 2,
    options = [];
elseif nargin == 3,
    options = varargin{1};
else 
    error('Incorrect number of input arguments.');
end
if ~isIQMmodel(model),
    error('Function only defined for IQMmodels.');
end
if iscell(data),
    for k=1:length(data),
        if ~isIQMmeasurement(data{k}),
            error('Error in the data input argument (IQMmeasurement required).');
        end
    end
elseif isIQMmeasurement(data),
    data = {data};
elseif isnumeric(data),
    time = data;
    if length(time) == 1,
        time = [0:time/1000:time];
    end
    data = {createDummyMeasurement(time,model)};
else
    error('Error in the data input argument (IQMmeasurement required).');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE DUMMY PROJECT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = IQMprojectSB();                          % create empty project
p = IQMupdatemodel(p,model);               % add model
e = IQMexperiment(); es = struct(e); es.name = 'Empty Experiment'; e = IQMexperiment(es);
p = IQMupdateexperiment(p,e); % add empty experiment
for k=1:length(data),
    p = IQMupdatemeasurement(p,1,data{k}); % add measurements
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine parameters changed by events
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters that are changed by events are not allowed to be changed by the user during tuning
ms = struct(model);
pnames = {ms.parameters.name};
% collect all assignment variables in all events that are parameters
apnames = {};
for k=1:length(ms.events),
    for k2=1:length(ms.events(k).assignment),
        vname = ms.events(k).assignment(k2).variable;
        if ~isempty(strmatchIQM(vname,pnames,'exact')),
            apnames{end+1} = vname;
        end
    end
end
apnames = unique(apnames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALL IQMmanualtuning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ptuned = modeltuningIQM(p,1,apnames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET TUNED MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modeltuned = IQMgetmodel(ptuned,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINSIHED => RETURN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a dummy measurement 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dummymeas] = createDummyMeasurement(time,model)
% Simulate the model over given time-vector
simdata = IQMPsimulate(model,time);
% create IQMmeasurement structure
ms = struct(IQMmeasurement);
ms.name = 'Nominal model';
ms.notes = '';
ms.time = time;
for k=1:length(simdata.states),
    ms.data(end+1).name = simdata.states{k};
    ms.data(end).notes = '';
    ms.data(end).values = simdata.statevalues(:,k);
    x = NaN(size(simdata.statevalues(:,k)));
    ms.data(end).maxvalues = x;
    ms.data(end).minvalues = x;
end
for k=1:length(simdata.variables),
    ms.data(end+1).name = simdata.variables{k};
    ms.data(end).notes = '';
    ms.data(end).values = simdata.variablevalues(:,k);
    x = NaN(size(simdata.variablevalues(:,k)));
    ms.data(end).maxvalues = x;
    ms.data(end).minvalues = x;
end
dummymeas = IQMmeasurement(ms);
return