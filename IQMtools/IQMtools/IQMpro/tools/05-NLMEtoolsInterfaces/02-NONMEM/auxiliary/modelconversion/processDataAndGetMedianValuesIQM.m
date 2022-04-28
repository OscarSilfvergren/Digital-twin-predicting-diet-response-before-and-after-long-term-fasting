function [ covariateMedianNames,covariateMedianValues,covariateCATNames,covariateCATValues,dataheader,dataCSV ] = processDataAndGetMedianValuesIQM( oldpath,dataRelPathFromProject,dataFileName,dataHeaderIdent,SILENT,COVcentering_covs,COVcentering_values  )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if datafile exists and csv file and load some information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataFile = fullfile(dataRelPathFromProject,dataFileName);
try
    dataheader = IQMloadCSVdataset(dataFile,1);
catch
    cd(oldpath);
    error('Please check if the data file "%s" exists.',dataFile)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if length of header identical to dataHeaderIdent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(explodePCIQM(dataHeaderIdent,',')) ~= length(dataheader),
    cd(oldpath);
    error('Please check: The data header identifiers do not have the same length as the number of columns in the dataset.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load dataset and get header information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataCSV     = IQMloadCSVdataset(fullfile(dataRelPathFromProject,dataFileName));
dataheader  = dataCSV.Properties.VariableNames;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the NLME dataset for the minimal required columns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IQMcheckNLMEdatasetHeader(dataCSV);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine medians for covariates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine index of COV columns and their names
terms                   = explodePCIQM(dataHeaderIdent);
ixCOVs                  = strmatchIQM('COV',terms,'exact');

covariateMedianValues = [];
covariateMedianNames = {};
if ~isempty(ixCOVs),
    dataheaderCOVs      = dataheader(ixCOVs);
    
    % Determine index of ID column and ID name
    terms               = explodePCIQM(dataHeaderIdent);
    ixID                = strmatchIQM('ID',terms,'exact');
    dataheaderID        = dataheader(ixID);
    
    % Get covariate values for each individual
    allID               = unique(dataCSV.(dataheaderID{1}));
    allCOVs             = NaN(length(allID),length(ixCOVs));
    for k=1:length(allID),
        datak           = dataCSV(dataCSV.(dataheaderID{1})==allID(k),ixCOVs);
        allCOVs(k,:)    = table2array(datak(1,:));
    end
    
    % Determine median
    covariateMedianValues   = median(allCOVs);
    covariateMedianNames    = dataheaderCOVs;
    
    % Handle custom centering values
    for k=1:length(COVcentering_covs),
        ix = strmatchIQM(COVcentering_covs{k},covariateMedianNames,'exact');
        covariateMedianValues(ix) = COVcentering_values(k);
    end
    
    if ~SILENT,
        disp(' ')
        disp('Analysis of dataset for covariates - determine the centering values  ')
        disp('These are the median values, if not defined differently by the user.')
        disp(' Results:');
        for k=1:length(covariateMedianValues),
            disp(sprintf('   median(%s) = %g',covariateMedianNames{k},covariateMedianValues(k)));
        end
        disp('These values will be used to center the continuous covariates')
        disp(' ')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determines the categorical covariates and the elements they can take
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine index of COV columns and their names
terms                   = explodePCIQM(dataHeaderIdent);
ixCATs                  = strmatchIQM('CAT',terms,'exact');

covariateCATValues = {};
covariateCATNames = {};
if ~isempty(ixCATs),
    dataheaderCATs      = dataheader(ixCATs);
    
    % Get unique elements for each cat covariate and names ...
    for k=1:length(ixCATs),
        covariateCATNames{k} = dataheaderCATs{k};
        covariateCATValues{k} = unique(dataCSV.(dataheaderCATs{k}));
        if sum(isnan(unique(dataCSV.(dataheaderCATs{k})))) > 0,
            error('The categorical covariate "%s" contains NaN => please impute before running the parameter estimation.',dataheaderCATs{k});
        end
    end
    
    % Print information about reference values
    if ~SILENT,
        disp(' ');
        disp('The following values for the categorical covariates are used as reference values:');
        for k=1:length(covariateCATNames),
            disp(sprintf('\t%s%s: %d',covariateCATNames{k},char(32*ones(1,cellmaxlengthIQM(covariateCATNames)-length(covariateCATNames{k})+5)),covariateCATValues{k}(1)));
        end
        disp(' ');
    end
end