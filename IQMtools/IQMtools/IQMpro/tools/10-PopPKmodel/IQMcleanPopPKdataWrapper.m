function [dataCleaned] = IQMcleanPopPKdataWrapper(data,DOSENAME,OBSNAME,covNames,catNames,removeSUBJECT,removeREC,outputPath,catImputationValues)
% This function is a wrapper for different cleaning functions for standard
% popPK datasets. Standard here means single compound and single PK
% readout.
% 
% The following "cleaning" functions are called in sequence:
% 
% 1) IQMcleanRemoveRecordsSUBJECTs - removes manually selected records and subjects
%
% 2) IQMselectDataEvents - keeps only selected dose and observation
%
% 3) IQMcleanRemovePlaceboSubjects - removes placebo subjects
% 
% 4) IQMcleanRemoveSubjectsNoObservations - removes subjects without PK observations
% 
% 5) IQMcleanRemoveZeroDoses - removes doses with 0 amount
%
% 6) IQMcleanImputeCovariates - imputes covariates - continuous to the median, categorical to 99 or to provided imputation values
%
% 7) Check if MDV=1 present for observations without a reason in IGNORE column.
%    If yes, then throw an error.
%
% 8) Set all non-empty IGNORE observation records to MDV=1
% 
% [SYNTAX]
% [dataCleaned] = IQMcleanPopPKdataWrapper(data,DOSENAME,OBSNAME,covNames,catNames)
% [dataCleaned] = IQMcleanPopPKdataWrapper(data,DOSENAME,OBSNAME,covNames,catNames,removeSUBJECT,removeREC)
% [dataCleaned] = IQMcleanPopPKdataWrapper(data,DOSENAME,OBSNAME,covNames,catNames,removeSUBJECT,removeREC,outputPath)
% [dataCleaned] = IQMcleanPopPKdataWrapper(data,DOSENAME,OBSNAME,covNames,catNames,removeSUBJECT,removeREC,outputPath,catImputationValues)
%
% [INPUT]
% data:           Dataset in task specific dataset format.
% DOSENAME:       String with the NAME of the dosing event to consider.
% OBSNAME:        String with the NAME of the observation event to consider.
% covNames:       Cell-array with names of continuous covariates
% catNames:       Cell-array with names of categorical covariates
% removeSUBJECT:  Cell-matrix with 2 columns. First column contains the
%                 USUBJID identifiers of the subjects to be removed
%                 from the dataset. The second column contains strings,
%                 which define the reason why this subject is removed.
%                 Default: {} - nothing removed.
% removeREC:      Cell-matrix with 2 columns. First column contains the indices
%                 of the records to be removed from the dataset. The second
%                 column contains strings, which define the reason why this
%                 record is removed. The indices correspond to the number
%                 of the row in the dataset (-1 due to the header). These
%                 indices are also displayed as IX# in the individual plots
%                 when using the IQMexplorePKdataWrapper function.
%                 Default: {} - nothing removed.
% outputPath:     Path where to store the cleaning log files
%                 Default: '../Output/DataCleaning/';
% catImputationValues: Vector with values for imputation of categorical
%                 covariates. If not provided "99" will be used for all
%                 imputed categorical covariates. Vector needs same order
%                 and number of elements as catNames. 
%
% [OUTPUT]
% dataCleaned:    Cleaned dataset 

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

% Handle variable input arguments
if nargin<5 || nargin>9,
    error('Incorrect number of input arguments.');
end
if nargin<6,
    removeSUBJECT = {};
end
if nargin<7,
    removeREC = {};
end
if nargin<8,
    outputPath = '../Output/DataCleaning/';
end
if nargin<9,
    catImputationValues = [];
end

% Handle cell thing
if ischar(covNames),
    covNames = {covNames};
end
if ischar(catNames),
    catNames = {catNames};
end

% Handle cell thingy
if ischar(OBSNAME),
    OBSNAME = {OBSNAME};
end
if length(OBSNAME) > 1,
    error('OBSNAME can only contain one NAME.');
end
OBSNAME = OBSNAME{1};

if ischar(DOSENAME),
    DOSENAME = {DOSENAME};
end
if length(DOSENAME) > 1,
    error('DOSENAME can only contain one NAME.');
end
DOSENAME = DOSENAME{1};

% Remove manually selected records and subjects 
% Records are not actually removed ... the IGNORE entries are set to the
% provided reason instead. And MDV is set to 1.
% It is checked that this is only done on observation records (EVID=0). If
% done on dose records an error is thrown.
% Subjects are removed completely from the dataset, as this involves also
% removing doses and to be not tool specific (e.g. MONOLIX does not know an
% IGNORE statement) we need to ask the user to 
% either accept removal completely or handle it outside this function.
filename = [outputPath '/01_Manually_Selected_Records_and_Subjects'];
datax = IQMcleanRemoveRecordsSUBJECTs(data,removeSUBJECT,removeREC,filename);

% Keep only selected dose and observation
datax = IQMselectDataEvents(datax,{DOSENAME OBSNAME}); 

% Remove placebo subjects
filename = [outputPath '/02_PLACEBO_subjects'];
datax = IQMcleanRemovePlaceboSubjects(datax,filename);

% Remove subjects without PK observations
filename = [outputPath '/03_Subjects_without_observations'];
datax = IQMcleanRemoveSubjectsNoObservations(datax,filename);

% Remove doses with 0 amount
filename = [outputPath '/04_Dose_Records_zero_amount'];
datax = IQMcleanRemoveZeroDoses(datax,filename);

% Covariate imputation
filename = [outputPath '/05_Covariate_Imputation'];
if isempty(catImputationValues),
    catImputationValues = 99*ones(1,length(catNames));
end
datax = IQMcleanImputeCovariates(datax,covNames,catNames,catImputationValues,filename);

% Check if MDV=1 present for observations without a reason in IGNORE column.
% If yes, then throw an error.
ix_MDV_IGNORE_test = find(datax.EVID==0 & datax.MDV==1 & strcmp(datax.IGNORE,''));
if ~isempty(ix_MDV_IGNORE_test),
    error('MDV=1 observations present that do not have an IGNORE reason.');
end

% Set all non-empty IGNORE observation records to MDV=1
ix_IGNORE = find(datax.EVID==0 & ~strcmp(datax.IGNORE,''));
datax.MDV(ix_IGNORE) = 1;

% Assign output
dataCleaned = datax;

% Check length of original data with cleaned dataset. When the compliance
% mode is on all removals should have been handled before (e.g. by using
% the IGNORE column) and adequate specification of the data set spec.
% Meaning that nothing should actually be really removed from the data.
SETUP_PATHS_TOOLS_IQMLITE  % Get compliance mode setting
if height(data) ~= height(dataCleaned),
    text = sprintf('The IQMcleanPopPKdataWrapper function removed some rows in the dataset.\n');
    text = sprintf('%sThis is perfectly normal and it is supposed to do that (see help text for this function).\n',text);
    text = sprintf('%sBut if you would like to keep track of removed records in clinical projects for compliance\n',text);
    text = sprintf('%sreasons, then this is not acceptable and removals of single records should be handled using\n',text);
    text = sprintf('%sthe IGNORE column. Removal of subjects should be handled prior to the call to IQMcleanPopPKdataWrapper.\n',text);
    text = sprintf('%sThis message appears as a warning is compliance mode is "OFF" and as an error if it is "ON".\n',text);
    if COMPLIANCE_OUTPUT_MODE == 1,
        error(text);
    else
        warning(text);
    end
end
