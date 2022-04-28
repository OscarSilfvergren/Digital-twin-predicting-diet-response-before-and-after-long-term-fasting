function [doseeventstruct] = dosing2doseeventIQM(dosing)
% dosing2doseeventIQM: This function takes an IQMdosing object as
% input and returns a structure in which the dosing events, happening, are
% sorted according to the times at which they happen. Additionally, the
% output structure contains information about the dosing amount and the
% names of the parameters in a model to change in order to apply the
% dosing. Also lag and the other dosing type dependend informations are
% collected. Lag information is only collected if lag is numeric ... if
% non-numeric then it is assumed to be defined in the model already as
% parameter of variable and then there is no need to get it from the dosing
% scheme.
% 
% USAGE:
% ======
% [doseeventstruct] = dosing2doseeventIQM(dosing) 
%
% dosing: IQMdosing object
%
% Output Arguments:
% =================
% doseeventstruct: output argument with the following structure:
%   .time                   time of dosing application
%   .input 
%   .input.index            index of dosing input in dosing object to be applied
%   .input.timeparname      name of time parameter for this input
%   .input.parameternames   cell-array with other parameter names to change the values for
%   .input.parametervalues  vector with new parameter values (same order as
%                           parameternames cell-array)

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isIQMdosing(dosing),
    error('Input argument is not an IQMdosing object.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get dosing structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ds = struct(dosing);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get a dosing matrix and dose/time parameter 
% names as added to a model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dosingmatrix = [];
inames = {}; 
for k=1:length(ds.inputs),
    dosingmatrix = [dosingmatrix; ds.inputs(k).time(:) ds.inputs(k).D(:) k*ones(length(ds.inputs(k).D),1)];
    inames{k} = ds.inputs(k).name;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check if dosing times unique
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if length(ds.inputs(k).time) ~= length(unique(ds.inputs(k).time)),
        error('Dosing times in dosing object are not unique.');
    end
end
% sort the dosing matrix according to the time vector
dosingmatrix = sortrows(dosingmatrix,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each input in dosing scheme get info about 
% Tlag and the parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputinfoadditional = [];
for k=1:length(ds.inputs),
    inputinfoadditional(k).input = ds.inputs(k).name; % unused but just to have it in there when looking at this structure
    inputinfoadditional(k).index = k; % unused but just to have it in there when looking at this structure
    try
        inputinfoadditional(k).paramname = ds.inputs(k).parameters.name;
        inputinfoadditional(k).param = ds.inputs(k).parameters.value;
    catch
        inputinfoadditional(k).paramname = {};
        inputinfoadditional(k).param = [];
    end
    inputinfoadditional(k).Tlag = ds.inputs(k).Tlag;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert dosing matrix to structure
% dosingevent.time
% dosingevent.input.index
% dosingevent.input.timeparname
% dosingevent.input.parameternames
% dosingevent.input.parametervalues
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputstruct = struct('index',{},'timeparname',{},'parameternames',{},'parametervalues',{});
dosingevent = struct('time',{},'input',inputstruct);
times = unique(dosingmatrix(:,1));
for k=1:length(times),
    % find dosings for current time
    di = find(dosingmatrix(:,1)==times(k));
    % add info to dosingevent structure
    dosingevent(k).time = times(k);
    for k2=1:length(di),
        index = dosingmatrix(di(k2),3);
        dosingevent(k).input(k2).index = index;
        dosingevent(k).input(k2).timeparname = ['Time' lower(inames{index})]; % Name definition used in mergemoddosIQM
        % other parameters
        parameternames = {};
        parameternames{end+1} = ['Dose' lower(inames{index})];
        if ~isempty(inputinfoadditional(index).Tlag),
            % Only add Tlag if it is numeric
            if isnumeric(inputinfoadditional(index).Tlag),
                parameternames{end+1} = ['Tlag' lower(inames{index})];
            end
        end
        
        if ~isempty(inputinfoadditional(index).paramname),
            parameternames{end+1} = [inputinfoadditional(index).paramname '' lower(inames{index})];
        end

        parametervalues = [];
        parametervalues(end+1) = dosingmatrix(di(k2),2);
        if ~isempty(inputinfoadditional(index).Tlag),
            % Only add Tlag if it is numeric
            if isnumeric(inputinfoadditional(index).Tlag),
                parametervalues(end+1) = inputinfoadditional(index).Tlag;
            end
        end

        if ~isempty(inputinfoadditional(index).param),
            timedoseinput = dosingevent(k).time;
            timeINPUT = ds.inputs(dosingevent(k).input(k2).index).time;
            paramvalueINPUT = ds.inputs(dosingevent(k).input(k2).index).parameters.value;
            pvalue = paramvalueINPUT(timeINPUT==timedoseinput);
            parametervalues(end+1) = pvalue;
        end
        
        dosingevent(k).input(k2).parameternames = parameternames;
        dosingevent(k).input(k2).parametervalues = parametervalues;
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Return the result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
doseeventstruct = dosingevent;
return