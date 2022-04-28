function IQMcreateMONOLIXproject(model,dosing,data,projectPath,varargin)
% Creates a Monolix/MLXTRAN project from an IQMmodel and an IQMdosing
% scheme. 
%
% This function will check if underscores are present in parameter and
% covariatenames. This is not allowed in IQM Tools for Monolix conversion.
% The reason is that Lixoft is constantly changing the output syntax of the
% parameter estimates and only by not allowing underscores it is possible
% to make sense out of the output across Monolix versions.
%
% IMPORTANT: continuous covariates are always log-transformed, independent
% of the distribution of the parameter on which they are added. 
%
% [SYNTAX]
% [] = IQMcreateMONOLIXproject(model,dosing,data,projectPath)
% [] = IQMcreateMONOLIXproject(model,dosing,data,projectPath,options)
% [] = IQMcreateMONOLIXproject(model,dosing,data,projectPath,options,parameterOrder)
%
% [INPUT]
% model:            IQMmodel - with additional annotation:
%                       Variables OUTPUT1-N where 1-N matches the YTYPE definitions in the dataset.
%                       <estimate> as comment on parameters to be estimated.
%                       <regression> as comment on parameters to be obtained from dataset.
% dosing:           IQMdosing object (or empty [] if no INPUT defined in model)
% data:             Structure with following fields:
%       data.dataRelPathFromProject:    path to data file - relative to the
%                                       projectPath folder.
%       data.dataFileName:              data file filename
%       data.dataHeaderIdent:           String with datafile header identifiers (example: 'ID,TIME,Y,MDV,EVID,AMT,TINF,ADM,YTYPE,COV,COV,CAT') 
% projectPath:      String with the path/foldername to which the project files are to be written (example: 'FIT_01' or 'Models/FITS/FIT_01') 
% parameterOrder:   Used to reorder parameters (used by the popPK workflow, do not use otherwise)
%
% options:      Structure with following fields (all optional with default settings):
%       options.POPestimate:            Vector with 0 and 1 entries. 1 if pop parameter is estimated, 0 if not. Default or []: => all are estimated
%       options.POPvalues0:             Vector with pop parameter initial values. Default or []: => values stored in model and dosing scheme
%       options.IIVdistribution:        Cell-array with information about parameter distribution. L (lognormal), N (normal), G (logit)
%                                       Example: {'L' 'L' 'L' 'L' 'N' 'L' 'L' 'L'}. Default or {}: => use lognormal for all
%       options.IIVestimate:            Vector with 0 and 1 entries. 1 if random effect is estimated, 0 if not. Default or []: => all are estimated
%                                       0: IIV not estimated (IIVvalues0 not used) 
%                                       1: IIV estimated (IIVvalues0 as starting guesses)
%                                       2: IIV not estimated but fixed on IIVvalues0 value
%       options.IIVvalues0:             Vector with random effect parameter
%                                       initial values. Default or []: => all set to 0.5
%                                       If IIV not estimated then defined initial guess not used but replaced by 0
%       options.errorModels:            String with definition of residual error models, comma separated for each output.
%                                       Possible values: const,prop,comb1. Example: 'comb1,prop', Default or '': => const for all outputs
%       options.errorParam0:            Vector allowing to pass initial guesses for error model parameters. Same order as error models. 
%                                       'const': a, 'prop': b, 'comb1': a,b
%       options.covarianceModel:        Definition of covariance model. String with cell-array text inside, grouping the parameters to consider having 
%                                       correlated random effects. Example: '{CL,Vc},{Q,Vp,KM}'. Default: 'diagonal'
%       options.covariateModel:         Definition of covariate model. Cell-array. Each element is a sub-cell-array. First element in sub-cell-array is the 
%                                       parameter to which to add the covariate, all following elements define the covariates as named in the dataset.
%                                       Example: '{CL,BMI0}, {Fsubcut,WT0}, {Vc,SEX,BMI0}'. Default: '' => no covariates
%                                       By default (and so far not changeable, the continuous covariates are all weighted by their median, determined from the dataset)
%                                       >>>Covariates can be added to all parameters for which not both IIVestimate and POPestimate are 0.
%       options.covariateModelValues:   Definition of covariate coefficients for the selected covariate model. 
%                                       Syntax is similar to options.covariateModel. It is a cell-array containing vectors instead of cell-arrays.
%                                       Each vector contains values for the covariate coefficients, matching the covariateModel definition order.
%                                       Example: if options.covariateModel = '{CL,BMI0,AGE0}, {Fsubcut,WT0}, {Vc,SEX,BMI0}'
%                                       Then: options.covariateModelValues = {[0.5,0], [0.75], [0,0]}
%                                       Defines the initial guesses for the covariate coefficients for BMI0 on CL to be 0.5, WT0 on Fsubcut to be 0.75, and the other ones are 0.
%                                       If not defined, all covariate coefficients start from 0 in the estimation.
%       options.COVestimate:            Same structure as options.covariateModelValues but with entries 0 or 1. 0 means not estimated, 1 means estimated.
%                                       By default all are estimated.
%                                       In the example above options.COVestimate = {[0,1], [1], [1,0]}   will estimate AE0 on CL, WT0 on Fsubcut, SEX on Vc.
%                                       The other coefficients will be kept fixed.
%       options.COVcentering.covs:      Cell-array with covariates that should be centered around a custom value. 
%       options.COVcentering.values:    Vector with centering values. 
%       options.Ntests:                 Doing robustness analysis - number of models to generate with different initial guesses (randomly generated based on POPvalues0)
%                                       Default: 1 (no robustness analysis, using initial guesses as provided)
%       options.std_noise_setting:      Standard deviation to use to add noise to the initial parameter guesses (default=0.5 (50%CV))
%                                       Normal:         Parameter_guess + std_noise_setting*Parameter_guess*randomNumbers(0-1)
%                                       Lognormal:      Parameter_guess * exp(std_noise_setting*randomNumbers(0-1))
%                                       Logitnormal:    Similar and between 0-1
%       options.algorithm.SEED:         Seed setting. Defualt: 123456
%       options.algorithm.K1:           First iterations. Default: 500
%       options.algorithm.K2:           Final iterations. Default: 200
%       options.algorithm.K1_AUTO:      Automatic first iteration number (0: off, 1: on). Default: 0
%       options.algorithm.K2_AUTO:      Automatic final iteration number (0: off, 1: on). Default: 0
%       options.algorithm.NRCHAINS:     Number of parallel chains. Default: 1
%       options.algorithm.LLsetting:    'linearization' (default) or 'importantsampling'
%       options.algorithm.FIMsetting:   'linearization' (default) or 'stochasticApproximation'
%       options.algorithm.INDIVparametersetting: 'conditionalMode' (default) ... others not considered for now. 
%       options.algorithm.startTime:    start time of integration. default: [] (not set).
%       options.SILENT:                 =0: do some output in the command window, =1: do no output in command window (default: 0)
%       options.keepProjectFolder:      =0: remover already existing folder, =1: keep it
%
% [OUTPUT]
% Generated MONOLIX Project at projectPath location
%
% [ASSUMPTIONS]
% - Always log transformed continuous covariates (centered by the median)
% - Always untransformed categorical covariates
% - FIM: always by computed by linearization
% - Individual parameters are always determined by conditional modes, 
% - Use Monolix default graphics settings

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

% Define Default Properties (Never changing)
projectName            = 'project';
resultsFolder          = 'RESULTS';

% Check input
if ischar(model),
    model = IQMmodel(model);
end
if ischar(dosing),
    dosing = IQMdosing(dosing);
end
if ~isIQMmodel(model),
    error('First input argument is not an IQMmodel.');
end
if ~isIQMdosing(dosing) && ~isempty(dosing),
    error('Second input argument is not an IQMdosing scheme.');
end
try
    dataRelPathFromProject = data.dataRelPathFromProject;
    dataFileName           = data.dataFileName;
    dataHeaderIdent        = data.dataHeaderIdent;
catch
    error('data input argument not defined correctly.');
end

% Removal of TIMEPOS in dataHeaderIdent: TIMEPOS only needed for NONMEM ...
dataHeaderIdent         = regexprep(dataHeaderIdent,'\<TIMEPOS\>','IGNORE');

% Handle variable input arguments
options = [];
parameterOrder = {};
if nargin==5,
    options = varargin{1};
    parameterOrder = {};
elseif nargin==6,
    options = varargin{1};
    parameterOrder = varargin{2};
end

% Handle optional arguments
try POPestimate                     = options.POPestimate;                      catch, POPestimate = [];                             end
try POPvalues0                      = options.POPvalues0;                       catch, POPvalues0 = [];                              end
try IIVdistribution                 = options.IIVdistribution;                  catch, IIVdistribution = {};                         end
try IIVestimate                     = options.IIVestimate;                      catch, IIVestimate = [];                             end
try IIVvalues0                      = options.IIVvalues0;                       catch, IIVvalues0 = [];                              end
try errorModels                     = options.errorModels;                      catch, errorModels = '';                             end
try errorParam0                     = options.errorParam0;                      catch, errorParam0 = [];                             end
try covarianceModel                 = options.covarianceModel;                  catch, covarianceModel = 'diagonal';                 end
try covariateModel                  = options.covariateModel;                   catch, covariateModel = '';                          end
try covariateModelValues            = options.covariateModelValues;             catch, covariateModelValues = {};                    end
try COVestimate                     = options.COVestimate;                      catch, COVestimate = {};                             end
try COVcentering_covs               = options.COVcentering.covs;                catch, COVcentering_covs = {};                       end
try COVcentering_values             = options.COVcentering.values;              catch, COVcentering_values = [];                     end
try SEED                            = options.algorithm.SEED;                   catch, SEED = 123456;                                end
try K1                              = options.algorithm.K1;                     catch, K1 = 500;                                     end
try K2                              = options.algorithm.K2;                     catch, K2 = 200;                                     end
try K1_AUTO                         = options.algorithm.K1_AUTO;                catch, K1_AUTO = 0;                                  end
try K2_AUTO                         = options.algorithm.K2_AUTO;                catch, K2_AUTO = 0;                                  end
try NRCHAINS                        = options.algorithm.NRCHAINS;               catch, NRCHAINS = 1;                                 end
try SILENT                          = options.SILENT;                           catch, SILENT = 0;                                   end
try INDIVparametersetting           = options.algorithm.INDIVparametersetting;  catch, INDIVparametersetting = 'conditionalMode';    end
try LLsetting                       = options.algorithm.LLsetting;              catch, LLsetting = 'linearization';                  end
try FIMsetting                      = options.algorithm.FIMsetting;             catch, FIMsetting = 'linearization';                 end
try startTime                       = options.algorithm.startTime;              catch, startTime = [];                               end
try keepProjectFolder               = options.keepProjectFolder;                catch, keepProjectFolder = 0;                        end   
try Ntests                          = options.Ntests;                           catch, Ntests = 1;                                   end
try std_noise_setting               = options.std_noise_setting;                catch, std_noise_setting = 0.5; options.std_noise_setting = 0.5;    end

% Handle cell requirements
if ~iscell(COVcentering_covs),
    COVcentering_covs = {COVcentering_covs};
end

options.covariateModel = covariateModel;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Handle first the case with robustness analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Ntests>1,
    % Check if POPvalues0 and POPestimate defined
    if isempty(POPestimate) || isempty(POPvalues0),
        error('When doing robustness analysis, please define options.POPvalues0 and options.POPestimate!');
    end
    
    % Define IIVdistribution if still empty (default: all 'L')
    if isempty(IIVdistribution),
        for k=1:length(POPvalues0),
            IIVdistribution{k} = 'L';
        end
    end
    
    % Sample Ntests new POPvalues0 for the ones that are estimated using
    % std_noise_setting as standard deviation
    
    % Allocating variable
    POPvalues0_sampled                              = POPvalues0(ones(1,Ntests),:);
    
    % Sampling normally distributed (IIV) parameters - which are also estimated on a population level
    ix_normal_sampled                               = find(strcmp(IIVdistribution,'N').*POPestimate);
    POPvalues0_sampled(:,ix_normal_sampled)         = POPvalues0(ones(1,Ntests),ix_normal_sampled) + std_noise_setting*POPvalues0(ones(1,Ntests),ix_normal_sampled).*randn(Ntests,length(ix_normal_sampled));

    % Sampling log normally distributed (IIV) parameters - which are also estimated on a population level
    ix_lognormal_sampled                            = find(strcmp(IIVdistribution,'L').*POPestimate);
    MU                                              = log(POPvalues0(ones(1,Ntests),ix_lognormal_sampled));
    XXX                                             = MU + std_noise_setting.*randn(Ntests,length(ix_lognormal_sampled));
    POPvalues0_sampled(:,ix_lognormal_sampled)      = exp(XXX);
    
    % Sampling logit normally distributed parameters - which are also estimated on a population level
    ix_logitnormal_sampled                          = find(strcmp(IIVdistribution,'G').*POPestimate);
    MU                                              = log(POPvalues0(ones(1,Ntests),ix_logitnormal_sampled)./(1-POPvalues0(ones(1,Ntests),ix_logitnormal_sampled)));
    XXX                                             = MU + std_noise_setting.*randn(Ntests,length(ix_logitnormal_sampled));
    POPvalues0_sampled(:,ix_logitnormal_sampled)    = exp(XXX)./(1+exp(XXX));
    
    % Clean folder
    try rmdir(projectPath,'s'); catch, end
    
    % Create Ntests different models in the projectPath/MODEL_01/02, ... folders
    for k=1:Ntests,
        % Setup new project creation stuff
        modelK                          = model;
        dosingK                         = dosing;
        dataK                           = data;
        dataK.dataRelPathFromProject    = ['../' data.dataRelPathFromProject];
        projectPathK                    = [projectPath sprintf('/MODEL_%s',preFillCharIQM(k,length(num2str(Ntests)),'0'))];
        optionsK                        = options;
        optionsK                        = rmfield(optionsK,'Ntests');
        optionsK                        = rmfield(optionsK,'std_noise_setting');
        optionsK.POPvalues0             = POPvalues0_sampled(k,:);
        IQMcreateMONOLIXproject(modelK,dosingK,dataK,projectPathK,optionsK)
    end
    
    % Ready, return
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert covariate model into different syntax
% '{CL,BMI0}, {Fsubcut,WT0}, {Vc,SEX,BMI0}'
% to
% {{'CL','BMI0'}, {'Fsubcut','WT0'}, {'Vc','SEX','BMI0'}}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(covariateModel),
    terms = explodePCIQM(covariateModel,',','{','}');
    y = {};
    for k=1:length(terms),
        x = strrep(strtrim(terms{k}),' ','');
        x = strrep(x,'{','{''');
        x = strrep(x,'}','''}');
        x = strrep(x,',',''',''');
        y{k} = eval(x);
    end
    covariateModel = y;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check covariateModelValues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(covariateModel) && ~isempty(covariateModelValues),
    error('If you define covariateModelValues, you also need to define the covariateModel.');
end

if isempty(covariateModelValues),
    % Determine default covariateModelValues
    covariateModelValues = {};
    for k=1:length(covariateModel),
        covariateModelValues{k} = zeros(1,length(covariateModel{k})-1);
    end
else
    % Check correct length of covariateModelValues elements
    if length(covariateModel) ~= length(covariateModelValues),
        error('Number of elements in covariateModel and covariateModelValues needs to match.');
    end
    for k=1:length(covariateModel),
        if length(covariateModel{k})-1 ~= length(covariateModelValues{k}),
            error('Length of single elements in covariateModel and covariateModelValues needs to match (covariateModelValues elements being one shorter).');
        end            
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check COVestimate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(covariateModel) && ~isempty(COVestimate),
    error('If you define COVestimate, you also need to define the covariateModel.');
end

if isempty(COVestimate),
    % Determine default COVestimate - all are estimates
    COVestimate = {};
    for k=1:length(covariateModel),
        COVestimate{k} = ones(1,length(covariateModel{k})-1);
    end
else
    % Check correct length of COVestimate elements
    if length(covariateModel) ~= length(COVestimate),
        error('Number of elements in covariateModel and COVestimate needs to match.');
    end
    for k=1:length(covariateModel),
        if length(covariateModel{k})-1 ~= length(COVestimate{k}),
            error('Length of single elements in covariateModel and COVestimate needs to match (COVestimate elements being one shorter).');
        end            
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create project and results folder
% Change into project path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning off
oldpath = pwd;
if ~keepProjectFolder,
    try, rmdir(projectPath,'s'); catch, end
end
mkdir(projectPath); cd(projectPath)
mkdir(resultsFolder);
warning on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load dataset and get column names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    dataCSV     = IQMloadCSVdataset(fullfile(dataRelPathFromProject,dataFileName));
catch
    error('Trouble loading the data file. Please check if data.dataRelPathFromProject has been defined correctly.');
end
dataheader  = dataCSV.Properties.VariableNames;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the NLME dataset for the minimal required columns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IQMcheckNLMEdatasetHeader(dataCSV);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine medians for covariates
% Also handle in case custom centering values are defined.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine index of COV columns and their names
terms   = explodePCIQM(dataHeaderIdent);
ixCOVs  = strmatchIQM('COV',terms,'exact');
if ~isempty(ixCOVs),
    covariateMedianNames            = dataheader(ixCOVs);

    % Get median values
    unique_covs                     = unique(dataCSV(:,{'ID' covariateMedianNames{:}}));
    covariateMedianValues           = nanmedianIQM(table2array(unique_covs(:,2:end)));

    % Handle custom centering values
    for k=1:length(COVcentering_covs),
        ix                          = strmatchIQM(COVcentering_covs{k},covariateMedianNames,'exact');
        covariateMedianValues(ix)   = COVcentering_values(k);
    end
    
    if ~SILENT, 
        disp(' ')
        disp('==================================================================');
        disp('Analysis of dataset for covariates - determine the centering values  ')
        disp('These are the median values, if not defined differently by the user.')
        disp(' Results:');
        for k=1:length(covariateMedianValues),
            disp(sprintf('   median(%s) = %g',covariateMedianNames{k},covariateMedianValues(k)));
        end
        disp('These values will be used to center the continuous covariates')
        disp('==================================================================');
        disp(' ')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get model structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ms = struct(model);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create MLXTRAN model
% Provide the correct regression parameter names as ordered in the dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
regressionParametersData = dataCSV.Properties.VariableNames(strmatchIQM('X',explodePCIQM(data.dataHeaderIdent),'exact'));
[modelFileName, modelInfo] = IQMcreateMLXTRANfile(model,dosing,'model',SILENT,regressionParametersData,startTime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If needed, reorder parameters
% Reordering of parameters happens only for parameters to be estimated and
% should never be used by a user. Only by the popPK workflow.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(parameterOrder),
    % First load the MLXTRAN file and change the order of the parameters
    ptext = sprintf('%s, ',parameterOrder{:});
    ptext = ptext(1:end-2);
    parameter_def_new = sprintf('parameter = {%s}\r',ptext);
    % Read file to be changed
    content = fileread(modelFileName);
    % Change
    content = regexprep(content,'(parameter = {[^\n]+)',parameter_def_new);
    % Save file again
    fid = fopen(modelFileName,'w');
    fprintf(fid,'%s',content);
    fclose(fid);
    % Reorder parameters in modelInfo.param_est
    ix_reorder = [];
    for k=1:length(parameterOrder),
        ix_reorder(k) = strmatchIQM(parameterOrder{k},{modelInfo.param_est.name},'exact');
    end
    modelInfo.param_est = modelInfo.param_est(ix_reorder);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Info text
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~SILENT,
    disp(' ')
    disp('==================================================================');
    disp('==================================================================');
    disp('== Start of creation of Monolix project.mlxtran file');
    disp('==================================================================');
    disp('==================================================================');
    disp(' ')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write out parameter names and values for information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~SILENT,
    disp('Parameters selected to be estimated and their values in the model (this order):');
    disp('===============================================================================');
    for k=1:length(modelInfo.param_est),
        fprintf('\t%d)\t%s%s: %g\n',k,modelInfo.param_est(k).name,char(32*ones(1,15-length(modelInfo.param_est(k).name))),modelInfo.param_est(k).value0(1));
    end
    disp('==================================================================');
    disp(' ');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if datafile exists and csv file and load some information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataFile = [dataRelPathFromProject '/' dataFileName];
try
    dataheader = IQMloadCSVdataset(dataFile,1);
catch
    cd(oldpath);
    error('Please check if the data file "%s" exists.',dataFile)
end
% Check if length of header identical to dataHeaderIdent
if length(explodePCIQM(dataHeaderIdent,',')) ~= length(dataheader),
    cd(oldpath);
    error('Please check: The data header identifiers do not have the same length as the number of columns in the dataset.')
end

% Determine continuous and categorical covariates
IDs = explodePCIQM(dataHeaderIdent,',');
covIDs = strmatchIQM('COV',upper(IDs));
covNames = dataheader(covIDs);
catIDs = strmatchIQM('CAT',upper(IDs));
catNames = dataheader(catIDs);

% % Check that regression parameters correctly defined
% This chcek is not needed anymore, since regression parameters taken from
% dataset according to user specs. And check that regression parameters are
% available in the model is done above in the write out of the regression
% parameters during creation of the MLXTRAN model file.
% nrREGSmodel = length(modelInfo.param_reg);
% nrREGSdata  = length(strmatchIQM('X',upper(IDs)));
% if nrREGSmodel ~= nrREGSdata,
%     cd(oldpath);
%     error('Different numbers of regression parameters in model and in dataset.');
% end

% This printout is not needed anymore, since regression parameters are not
% taken from the model but from the dataset and the header type information
% when creating a MONOLIX project.
% % Print table with regression parameters model and data
% REGAnamesData  = dataheader(strmatchIQM('X',upper(IDs)));
% if ~SILENT,
%     disp(' ');
%     disp('Please check that the following matches regression parameters in data and model do make sense:');
%     for k=1:length(REGAnamesData),
%         fprintf('\t%s%s: %s\n',REGAnamesData{k},char(32*ones(1,8-length(REGAnamesData{k}))),modelInfo.param_reg(k).name)
%     end
% end
% % Print table with header names and identifiers
% IDs = explodePCIQM(dataHeaderIdent,',');

if ~SILENT,
    disp(' ');
    disp('Please check that the following matches between data header and identifiers do make sense:');
    for k=1:length(dataheader),
        fprintf('\t%s%s: %s\n',dataheader{k},char(32*ones(1,8-length(dataheader{k}))),IDs{k})
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check parameter and covariate names ... '_' not allowed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parNamesCheck = [{modelInfo.param_est.name} covNames catNames];
text_error = '';
for k=1:length(parNamesCheck),
    if ~isempty(strfind(parNamesCheck{k},'_')),
        text_error = sprintf('%sUnderscores "_" are not allowed in parameter and covariate names: "%s".\n',text_error,parNamesCheck{k});
    end
end
if ~isempty(text_error),
    error(text_error);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check POPestimate thing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(POPestimate),
    POPestimate = ones(1,length(modelInfo.param_est));
end
if length(modelInfo.param_est) ~= length(POPestimate),
    cd(oldpath);
    error('Please make sure POPestimate is of same length as number of parameters to be estimated.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check POPvalues0 thing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(POPvalues0),
    POPvalues0 = [];
    for k=1:length(modelInfo.param_est),
        POPvalues0(k) = modelInfo.param_est(k).value0(1);
    end
end
if length(modelInfo.param_est) ~= length(POPvalues0),
    cd(oldpath);
    error('Please make sure POPvalues0 is of same length as number of parameters to be estimated.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check IIV distribution things
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(IIVdistribution),
    IIVdistribution = {};
    for k=1:length(modelInfo.param_est),
        IIVdistribution{k} = 'L';
    end
end

% Check contents
test = IIVdistribution;
for k=1:length(IIVdistribution),
    if ~ismember(test{k},{'L','N','G'}),
        cd(oldpath);
        error('Please make sure that only "N", "L", or "G" appear in the "IIVdistribution" variable.');
    end
end

% Check length
if length(IIVdistribution) ~= length(modelInfo.param_est),
    cd(oldpath);
    error('Please make sure that an equal number of IIVdistribution is defined as estimated parameters in the model.');
end

% Print table parameter names and IIV distributions
if ~SILENT,
    disp(' ');
    disp('Please check that the following matches between parameters and used IIV distributions are correct:');
    for k=1:length(modelInfo.param_est),
        if IIVdistribution{k} == 'L', dtext = 'logNormal'; end
        if IIVdistribution{k} == 'N', dtext = 'Normal'; end
        if IIVdistribution{k} == 'G', dtext = 'logitNormal'; end
        fprintf('\t%s%s: %s\n',modelInfo.param_est(k).name,char(32*ones(1,15-length(modelInfo.param_est(k).name))),dtext)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check IIV estimation things
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(IIVestimate),
    IIVestimate = ones(1,length(modelInfo.param_est));
end
% Check length
if length(IIVestimate) ~= length(modelInfo.param_est),
    cd(oldpath);
    error('Please make sure that an equal number of IIVestimate is defined as estimated parameters in the model.');
end
% Print table parameter names and IIV esimations
if ~SILENT,
    disp(' ');
    disp('Please check that the following matches between parameters and used estimated IIVs are correct:');
    for k=1:length(modelInfo.param_est),
        if IIVestimate(k) == 0,
            fprintf('\t%s%s: IIV NOT ESTIMATED (kept on 0)\n',modelInfo.param_est(k).name,char(32*ones(1,15-length(modelInfo.param_est(k).name))));
        elseif IIVestimate(k) == 1,
            fprintf('\t%s%s: IIV ESTIMATED\n',modelInfo.param_est(k).name,char(32*ones(1,15-length(modelInfo.param_est(k).name))));
        elseif IIVestimate(k) == 2,
            fprintf('\t%s%s: IIV NOT ESTIMATED (kept on initial value)\n',modelInfo.param_est(k).name,char(32*ones(1,15-length(modelInfo.param_est(k).name))));
        end
    end
    disp(' ');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check IIVvalues0 thing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(IIVvalues0),
    IIVvalues0 = 0.5*ones(1,length(modelInfo.param_est));
end
if length(modelInfo.param_est) ~= length(IIVvalues0),
    cd(oldpath);
    error('Please make sure IIVvalues0 is of same length as number of parameters to be estimated.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check residual error things
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(errorModels),
    errorModels = '';
    for k=1:length(modelInfo.outputs),
        errorModels = sprintf('%sconst,',errorModels);
    end
    errorModels = errorModels(1:end-1);
end
test = errorModels;
test = strtrim(strrep(strrep(strrep(strrep(strrep(test,'exp',''),'const',''),'prop',''),'comb1',''),',',''));
if ~isempty(test),
    cd(oldpath);
    error('Please make sure that only "const", "prop", "comb1", or "exp" appear in the "errorModels" variable.');
end
% Check length
errors = explodePCIQM(errorModels,',');
if length(errors) ~= length(modelInfo.outputs),
    cd(oldpath);
    error('Please make sure that an equal number of errorModels is defined as outputs in the model.');
end
% Print table parameter names and IIV distributions
if ~SILENT,
    disp(' ');
    disp('Please check that the following matches between outputs and used residual error models are correct:');
    for k=1:length(modelInfo.outputs),
        fprintf('\t%s%s: %s\n',modelInfo.outputs(k).formula,char(32*ones(1,15-length(modelInfo.outputs(k).formula))),errors{k})
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Handle empty errorParam0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(errorParam0),
    terms = explodePCIQM(errorModels);
    for k=1:length(terms),
        if strcmp(lower(terms{k}),'const'),
            errorParam0(end+1) = 1;
        elseif strcmp(lower(terms{k}),'prop'),
            errorParam0(end+1) = 0.3;
        elseif strcmp(lower(terms{k}),'comb1'),
            errorParam0(end+1) = 1;
            errorParam0(end+1) = 0.3;          
        elseif strcmp(lower(terms{k}),'exp'),
            errorParam0(end+1) = 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check errorParam0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
terms = explodePCIQM(errorModels);
nrneededelements = 0;
for k=1:length(terms),
    if strcmpi(terms{k},'const'),
        nrneededelements = nrneededelements+1;
    elseif strcmpi(terms{k},'prop'),
        nrneededelements = nrneededelements+1;
    elseif strcmpi(terms{k},'comb1'),
        nrneededelements = nrneededelements+2;
    elseif strcmpi(terms{k},'exp'),
        nrneededelements = nrneededelements+1;
    end
end
if length(errorParam0) ~= nrneededelements,
    error('Incorrect number of elements in options.errorParam0.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check covariance model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(covarianceModel),
    covarianceModel = 'diagonal';
elseif ~strcmp(covarianceModel,'diagonal'),
    % Need to check that none of the parameters for which no IIV is estimated is used in the covarianceModel
    param_est_noIIV = {modelInfo.param_est(find(~IIVestimate)).name};
    for k=1:length(param_est_noIIV),
        if ~isempty(regexp(covarianceModel,['\<' param_est_noIIV{k} '\>'])),
            cd(oldpath);
            error('Please make sure none of the parameters for which no IIV is estimated/fixed to non zero is used in the covarianceModel settings.');
        end
    end
    % Check that all parameters in the covariance model actually are model parameters
    param = {modelInfo.param_est.name};
    test  = covarianceModel;
    for k=1:length(param),
        test = regexprep(test,['\<' param{k} '\>'],'');
    end
    test = strrep(test,'{','');
    test = strrep(test,'}','');
    test = strrep(test,',','');
    test = strtrim(test);
    if ~isempty(test),
        cd(oldpath);
        error('Please make sure that covarianceModel only contains parameter names that are set to <estimate> in the model.');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check LL setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(LLsetting),
    LLsetting = 'both';
end
if isempty(strmatchIQM(LLsetting,{'linearization','importantSampling','both'})),
    cd(oldpath);
    error('Please make sure LLsetting has one of the following values: "linearization", "importantSampling", "both"=""');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check FIM setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(FIMsetting),
    FIMsetting = 'linearization';
end
if isempty(strmatchIQM(FIMsetting,{'linearization','stochasticApproximation'})),
    cd(oldpath);
    error('Please make sure FIMsetting has one of the following values: "linearization" or "stochasticApproximation"');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check covariate model things
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First check that all first elements are estimated parameters of the model
param_est = {modelInfo.param_est.name};
for k=1:length(covariateModel),
    param = covariateModel{k}{1};
    if isempty(strmatchIQM(param,param_est,'exact')),
        cd(oldpath);
        error('Please make sure that all parameters for which covariates are defined are defined by <estimate> in the model.');
    end
end
% Second check that all defined covariates actually are covariates
covcatNames = [covNames catNames];
for k=1:length(covariateModel),
    for k2=2:length(covariateModel{k}),
        cov = covariateModel{k}{k2};
        if isempty(strmatchIQM(cov,covcatNames,'exact')),
            cd(oldpath);
            error(sprintf('Please make sure that all covariates, defined in covariateModel, are defined in the dataset\n   This error might be due to a categorical covariate having only single category.'));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPEN File
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen([projectName '.mlxtran'],'w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INFO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'; MLXTRAN PROJECT, created using IQM Tools\r\n');
fprintf(fid,'; Date: %s\r\n',datestr(now,'yyyy-mmm-DD HH:MM'));
fprintf(fid,'; By:   %s\r\n',usernameIQM());
fprintf(fid,'\r\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Placeholder for project information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'; ==PROJECT HEADER START===================================================\r\n');
fprintf(fid,'PROJECT_HEADER_PLACEHOLDER\r\n');
fprintf(fid,'; ==PROJECT HEADER END=====================================================\r\n');
fprintf(fid,'\r\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'; =============================================\r\n');
fprintf(fid,'DESCRIPTION:\r\n');
fprintf(fid,'; =============================================\r\n');
fprintf(fid,'\t%s\r\n',ms.name);
fprintf(fid,'\r\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'; =============================================\r\n');
fprintf(fid,'DATA:\r\n');
fprintf(fid,'; =============================================\r\n');
fprintf(fid,'\tpath = "%%MLXPROJECT%%/%s/",\r\n',dataRelPathFromProject);
fprintf(fid,'\tfile  ="%s",\r\n',dataFileName);
fprintf(fid,'\theaders = {%s},\r\n',dataHeaderIdent);
fprintf(fid,'\tcolumnDelimiter = ","\r\n');
fprintf(fid,'\r\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'; =============================================\r\n');
fprintf(fid,'VARIABLES:\r\n');
fprintf(fid,'; =============================================\r\n');
% Assume all covariates defined in dataset are used ...
% Continuous are always log transformed and centered by median
% Categorical are always kept as they are
text = '';
% Write out continuous first
for k=1:length(covNames),
    text = sprintf('%s\t%s,\r\n',text,covNames{k});
    covname = covNames{k};
    % Scale covariate by median value from dataset
    ixmedian = strmatchIQM(covname,covariateMedianNames,'exact');
    covname_weighted = sprintf('%s/%g',covname,covariateMedianValues(ixmedian));
    text = sprintf('%s\tt_%s = log(%s) [use=cov],\r\n',text,covname,covname_weighted);
end
% Write out categorical
for k=1:length(catNames),
    text = sprintf('%s\t%s [use=cov, type=cat],\r\n',text,catNames{k});
end
% Remove last comma and write text to file
fprintf(fid,text(1:end-3));
fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

% Create BETACOVNAMES and BETACOVTRANS information for header
BETACOVNAMES = {};
BETACOVTRANS = {};
for k=1:length(covariateModel),
    for k2=2:length(covariateModel{k}),
        if ~isempty(strmatchIQM(covariateModel{k}{k2},covNames,'exact')),
            BETACOVNAMES{end+1} = sprintf('beta_%s(%s)',covariateModel{k}{1},covariateModel{k}{k2});
            ixmedian = strmatchIQM(covariateModel{k}{k2},covariateMedianNames,'exact');
            BETACOVTRANS{end+1} = sprintf('log(cov/%g)',covariateMedianValues(ixmedian));
        end
    end
end

% Create BETACATNAMES and BETACATREFERENCE information for header
BETACATNAMES        = {};
BETACATREFERENCE    = [];
for k=1:length(covariateModel),
    for k2=2:length(covariateModel{k}),
        if ~isempty(strmatchIQM(covariateModel{k}{k2},catNames,'exact')),
            BETACATNAMES{end+1} = sprintf('beta_%s(%s)',covariateModel{k}{1},covariateModel{k}{k2});
            BETACATREFERENCE(end+1) = min(unique(dataCSV.(covariateModel{k}{k2})));
        end
    end
end

% Determine all categories for categorical covariates and store them as
% metadata in the header of the project.mlxtran file
CAT_CATEGORIES = '';
for k=1:length(catNames),
    x = unique(dataCSV.(catNames{k}));
    x = sprintf('%g,',x);
    x = ['[' x(1:end-1) ']'];
    CAT_CATEGORIES{k} = x;
end
if ~isempty(CAT_CATEGORIES),
    CAT_CATEGORIES = sprintf('%s,',CAT_CATEGORIES{:});
    CAT_CATEGORIES = [CAT_CATEGORIES(1:end-1)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INDIVIDUAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'; =============================================\r\n');
fprintf(fid,'INDIVIDUAL:\r\n');
fprintf(fid,'; =============================================\r\n');
% Write out parameters to estimate. Assume iiv=yes on all of them by default. 
% If no IIV desired then rather fix omega to 0.01.
text = '';
PARAM_TRANSNAME_STRING = {};
PARAM_INVTRANSNAME_STRING = {};

for k=1:length(modelInfo.param_est),
    if IIVdistribution{k} == 'L', 
        dtext = 'logNormal'; 
        PARAM_INVTRANSNAME_STRING{k} = 'log(psi)';
        PARAM_TRANSNAME_STRING{k} = 'exp(phi)';        
    end
    if IIVdistribution{k} == 'N', 
        dtext = 'Normal';
        PARAM_INVTRANSNAME_STRING{k} = '(psi)';
        PARAM_TRANSNAME_STRING{k} = '(phi)';
    end
    if IIVdistribution{k} == 'G', 
        dtext = 'logitNormal'; 
        PARAM_INVTRANSNAME_STRING{k} = 'log(psi./(1-psi))';
        PARAM_TRANSNAME_STRING{k} = 'exp(phi)./(1+exp(phi))';        
    end
    % check if IIV 
    if IIVestimate(k) == 0,
        iiv='no';
    else
        iiv='yes';
    end
    % check for covariates to use
    param = modelInfo.param_est(k).name;
    covs = {};
    for k2=1:length(covariateModel),
        if strcmp(param,covariateModel{k2}{1}),
            covs = covariateModel{k2}(2:end);
        end
    end
    % Attach "t_" to continuous covariate names, keep categorical covariate names same
    for k2=1:length(covs),
        if ~isempty(strmatchIQM(covs{k2},covNames,'exact')),
            covs{k2} = ['t_' covs{k2}];
        end
    end
    % Write it out
    if isempty(covs),
        text = sprintf('%s\t%s = {distribution=%s, iiv=%s},\r\n',text,param,dtext,iiv);
    else
        % Create cov text
        covText = '';
        for k2=1:length(covs),
            covText = sprintf('%s%s,',covText,covs{k2});
        end
        covText = covText(1:end-1);
        text = sprintf('%s\t%s = {distribution=%s, covariate={%s}, iiv=%s},\r\n',text,param,dtext,covText,iiv);
        

    end        
end
% Remove last comma and write text to file
fprintf(fid,text(1:end-3));
fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CORRELATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp(covarianceModel,'diagonal'),
    fprintf(fid,'; =============================================\r\n');
    fprintf(fid,'CORRELATION:\r\n');
    fprintf(fid,'; =============================================\r\n');
    fprintf(fid,'\tcorrelationIIV = {%s}\r\n',covarianceModel);
    fprintf(fid,'\r\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STRUCTURAL_MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'; =============================================\r\n');
fprintf(fid,'STRUCTURAL_MODEL:\r\n');
fprintf(fid,'; =============================================\r\n');
fprintf(fid,'\tfile = "mlxt:%s",\r\n',strrep(modelFileName,'.txt',''));
fprintf(fid,'\tpath = "%%MLXPROJECT%%",\r\n');
fprintf(fid,'\toutput = {');
for k=1:length(modelInfo.outputs)-1,
    fprintf(fid,'%s, ',modelInfo.outputs(k).formula);
end
fprintf(fid,'%s}',modelInfo.outputs(end).formula);
fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBSERVATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'; =============================================\r\n');
fprintf(fid,'OBSERVATIONS:\r\n');
fprintf(fid,'; =============================================\r\n');
% Only consider "continuous" observations with IQM conversion 
errors = explodePCIQM(errorModels,',');
text = '';
for k=1:length(modelInfo.outputs),
    if strcmp(errors{k},'const'), errorModel = 'constant'; end
    if strcmp(errors{k},'prop'), errorModel = 'proportional'; end
    if strcmp(errors{k},'comb1'), errorModel = 'combined1'; end
    if strcmp(errors{k},'exp'), errorModel = 'exponential'; end
    text = sprintf('%s\ty%d = {type=continuous, prediction=%s, error=%s},\r\n',text,k,modelInfo.outputs(k).formula,errorModel);
end
% Remove last comma and write text to file
fprintf(fid,text(1:end-3));
fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TASKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'; =============================================\r\n');
fprintf(fid,'TASKS:\r\n');
fprintf(fid,'; =============================================\r\n');
fprintf(fid,'\t; settings\r\n');
fprintf(fid,'\tglobalSettings={\r\n');
fprintf(fid,'\t\twithVariance=no,\r\n'); % Always estimate standard deviations of IIVs
% Use default graphics settings
fprintf(fid,'\t\tsettingsGraphics="%%MLXPROJECT%%/project_graphics.xmlx",\r\n');
fprintf(fid,'\t\tsettingsAlgorithms="%%MLXPROJECT%%/project_algorithms.xmlx",\r\n');
fprintf(fid,'\t\tresultFolder="%%MLXPROJECT%%/%s"},\r\n',resultsFolder);
fprintf(fid,'\t; workflow\r\n');
fprintf(fid,'\testimatePopulationParameters(\r\n');
fprintf(fid,'\t\tinitialValues={\r\n');
% write out population parameter initial values
for k=1:length(modelInfo.param_est),
    method = '';
    if POPestimate(k) == 0,
        method = '[method=FIXED]';
    end
    fprintf(fid,'\t\t\tpop_{%s} = %g %s,\r\n',modelInfo.param_est(k).name,POPvalues0(k),method);
end

% write out covariate coefficient initial guesses
for k1=1:length(covariateModel),
    for k2=2:length(covariateModel{k1}),
        covarvalue = covariateModelValues{k1}(k2-1);
        if COVestimate{k1}(k2-1),
            method = '';
        else
            method = '[method=FIXED]';
        end
        ix = strmatchIQM(covariateModel{k1}{k2},covNames,'exact');
        if isempty(ix),
            fprintf(fid,'\t\t\tbeta_{%s,%s} = %g %s,\r\n',covariateModel{k1}{1},covariateModel{k1}{k2},covarvalue,method);
        else
            fprintf(fid,'\t\t\tbeta_{%s,t_%s} = %g %s,\r\n',covariateModel{k1}{1},covariateModel{k1}{k2},covarvalue,method);
        end
    end
end

% write out residual error model
errors = explodePCIQM(errorModels,',');
count = 1;
for k=1:length(errors),
    if strcmp(errors{k},'const'), 
        fprintf(fid,'\t\t\ta_y%d = %g,\r\n',k,errorParam0(count));
        count = count + 1;
    end
    if strcmp(errors{k},'prop'), 
        fprintf(fid,'\t\t\tb_y%d = %g,\r\n',k,errorParam0(count));
        count = count + 1;
    end
    if strcmp(errors{k},'comb1'), 
        fprintf(fid,'\t\t\ta_y%d = %g,\r\n',k,errorParam0(count));
        count = count + 1;
        fprintf(fid,'\t\t\tb_y%d = %g,\r\n',k,errorParam0(count));
        count = count + 1;
    end
    if strcmp(errors{k},'exp'), 
        fprintf(fid,'\t\t\ta_y%d = %g,\r\n',k,errorParam0(count));
        count = count + 1;
    end
end
% write out population parameter initial values
text = '';
for k=1:length(modelInfo.param_est),
    if IIVestimate(k)==1,
        value0 = IIVvalues0(k);
        text = sprintf('%s\t\t\tomega_{%s} = %g,\r\n',text,modelInfo.param_est(k).name,value0);
    elseif IIVestimate(k)==2,
        value0 = IIVvalues0(k);
        text = sprintf('%s\t\t\tomega_{%s} = %g [method=FIXED],\r\n',text,modelInfo.param_est(k).name,value0);
    end
end
fprintf(fid,text(1:end-3));
fprintf(fid,'\r\n');

fprintf(fid,'\t\t} ),\r\n');
if strcmp(FIMsetting,'linearization'),
    fprintf(fid,'\testimateFisherInformationMatrix( method={linearization} ),\r\n');
else
    fprintf(fid,'\testimateFisherInformationMatrix( method={stochasticApproximation} ),\r\n');
end    
fprintf(fid,'\testimateIndividualParameters( method={%s} ),\r\n',INDIVparametersetting);
if strcmp(LLsetting,'linearization'),
    fprintf(fid,'\testimateLogLikelihood(method={linearization}),\r\n');
elseif strcmp(LLsetting,'importantSampling'),
    fprintf(fid,'\testimateLogLikelihood(method={importantSampling}),\r\n');
elseif strcmp(LLsetting,'both'),
    fprintf(fid,'\testimateLogLikelihood(method={importantSampling,linearization}),\r\n');
end
fprintf(fid,'\tdisplayGraphics()');
fprintf(fid,'\r\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLOSE File
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Project Header with Metadata
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PROJECT_INFO_TEXT = '';

% Data location
DATA_info = sprintf('; DATA                = ''%s''\r\n',strrep(fullfile(dataRelPathFromProject,dataFileName),'\','/'));
PROJECT_INFO_TEXT = sprintf('%s%s',PROJECT_INFO_TEXT,DATA_info);

% DOSINGTYPES
ds = struct(dosing);
DOSINGTYPES = {ds.inputs.type};
x = sprintf('%s,',DOSINGTYPES{:});
DOSINGTYPES_info = sprintf('; DOSINGTYPES         = ''%s''\r\n',x(1:end-1));
PROJECT_INFO_TEXT = sprintf('%s%s',PROJECT_INFO_TEXT,DOSINGTYPES_info);

% covNames
x = sprintf('%s,',covNames{:});
COVNAMES_info = sprintf('; COVNAMES            = ''%s''\r\n',x(1:end-1));
PROJECT_INFO_TEXT = sprintf('%s%s',PROJECT_INFO_TEXT,COVNAMES_info);

% catNames
x = sprintf('%s,',catNames{:});
CATNAMES_info = sprintf('; CATNAMES            = ''%s''\r\n',x(1:end-1));
PROJECT_INFO_TEXT = sprintf('%s%s',PROJECT_INFO_TEXT,CATNAMES_info);

% CATCATEGORIES
CATCATEGORIES_info = sprintf('; CATCATEGORIES       = ''%s''\r\n',CAT_CATEGORIES);
PROJECT_INFO_TEXT = sprintf('%s%s',PROJECT_INFO_TEXT,CATCATEGORIES_info);

% Regression parameters
x = sprintf('%s,',regressionParametersData{:});
REGRESSNAMES_info = sprintf('; REGRESSIONNAMES     = ''%s''\r\n',x(1:end-1));
PROJECT_INFO_TEXT = sprintf('%s%s',PROJECT_INFO_TEXT,REGRESSNAMES_info);

% Outputs
x = cell(1,length(modelInfo.outputs));
for k=1:length(modelInfo.outputs),
    on = eval(strrep(modelInfo.outputs(k).name,'OUTPUT',''));
    x{on} = modelInfo.outputs(k).formula;
end
y = '';
for k=1:length(x),
    y = sprintf('%s%s,',y,x{k});
end
y = y(1:end-1);
OUTPUTS_info = sprintf('; OUTPUTS             = ''%s''\r\n',y);
PROJECT_INFO_TEXT = sprintf('%s%s',PROJECT_INFO_TEXT,OUTPUTS_info);

% Error models
ERRORMODELS_info = sprintf('; ERRORMODELS         = ''%s''\r\n',errorModels);
PROJECT_INFO_TEXT = sprintf('%s%s',PROJECT_INFO_TEXT,ERRORMODELS_info);

% PARAMNAMES
x = sprintf('%s,',modelInfo.param_est.name);
PARAMNAMES_info = sprintf('; PARAMNAMES          = ''%s''\r\n',x(1:end-1));
PROJECT_INFO_TEXT = sprintf('%s%s',PROJECT_INFO_TEXT,PARAMNAMES_info);

% PARAMTRANS
x = sprintf('%s,',PARAM_TRANSNAME_STRING{:});
PARAMTRANS_info = sprintf('; PARAMTRANS          = ''%s''\r\n',x(1:end-1));
PROJECT_INFO_TEXT = sprintf('%s%s',PROJECT_INFO_TEXT,PARAMTRANS_info);

% PARAMINVTRANS
x = sprintf('%s,',PARAM_INVTRANSNAME_STRING{:});
PARAMINVTRANS_info = sprintf('; PARAMINVTRANS       = ''%s''\r\n',x(1:end-1));
PROJECT_INFO_TEXT = sprintf('%s%s',PROJECT_INFO_TEXT,PARAMINVTRANS_info);

% COVARIATENAMES
COVARIATENAMES = [covNames,catNames];
x = sprintf('%s,',COVARIATENAMES{:});
COVARIATENAMES_info = sprintf('; COVARIATENAMES      = ''%s''\r\n',x(1:end-1));
PROJECT_INFO_TEXT = sprintf('%s%s',PROJECT_INFO_TEXT,COVARIATENAMES_info);

% COVARIATESUSED
COVARIATESUSED = setdiff(explodePCIQM(strrep(strrep(options.covariateModel,'{',''),'}','')),param_est);
x = sprintf('%s,',COVARIATESUSED{:});
COVARIATESUSED_info = sprintf('; COVARIATESUSED      = ''%s''\r\n',x(1:end-1));
PROJECT_INFO_TEXT = sprintf('%s%s',PROJECT_INFO_TEXT,COVARIATESUSED_info);

% BETACOVNAMES
x = sprintf('%s,',BETACOVNAMES{:});
BETACOVNAMES_info = sprintf('; BETACOVNAMES        = ''%s''\r\n',x(1:end-1));
PROJECT_INFO_TEXT = sprintf('%s%s',PROJECT_INFO_TEXT,BETACOVNAMES_info);

% BETACOVTRANS
x = sprintf('%s,',BETACOVTRANS{:});
BETACOVTRANS_info = sprintf('; BETACOVTRANS        = ''%s''\r\n',x(1:end-1));
PROJECT_INFO_TEXT = sprintf('%s%s',PROJECT_INFO_TEXT,BETACOVTRANS_info);

% BETACATNAMES
x = sprintf('%s,',BETACATNAMES{:});
BETACATNAMES_info = sprintf('; BETACATNAMES        = ''%s''\r\n',x(1:end-1));
PROJECT_INFO_TEXT = sprintf('%s%s',PROJECT_INFO_TEXT,BETACATNAMES_info);

% BETACATREFERENCE
x = ''; for k=1:length(BETACATREFERENCE), x=sprintf('%s%g,',x,BETACATREFERENCE(k)); end
BETACATREFERENCE_info = sprintf('; BETACATREFERENCE    = ''%s''\r\n',x(1:end-1));
PROJECT_INFO_TEXT = sprintf('%s%s',PROJECT_INFO_TEXT,BETACATREFERENCE_info);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Replace PROJECT_HEADER_PLACEHOLDER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
content = fileread('project.mlxtran');
content = strrep(content,'PROJECT_HEADER_PLACEHOLDER',strtrim(PROJECT_INFO_TEXT));
fid = fopen('project.mlxtran','w');
fprintf(fid,'%s',content);
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do generate the default graphics settings so that 
% predictions.txt file is generated and included NPDE and meanPWRES
% Trick is to load project file and to add things and then to save the file
% again.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

copyfile(which('template_project_graphics.xmlx'),'project_graphics.xmlx')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate project_algorithms.xmlx file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% only writing out K1, K2, Number of Chains, Seed, and auto settings for K1,K2,Number of chains
fid = fopen('project_algorithms.xmlx','w');
fprintf(fid,'<monolix>\n');
fprintf(fid,'	<algorithms seed="%d">\n',SEED);
fprintf(fid,'		<populationParameters>\n');
fprintf(fid,'			<vna value="%d,%d"/>\n',K1,K2);
fprintf(fid,'			<iop_Kauto value="%d,%d"/>\n',K1_AUTO,K2_AUTO);
fprintf(fid,'			<nmc value="%d"/>\n',NRCHAINS);
fprintf(fid,'		</populationParameters>\n');
fprintf(fid,'	</algorithms>\n');
fprintf(fid,'</monolix>\n');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change out of project path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(oldpath);