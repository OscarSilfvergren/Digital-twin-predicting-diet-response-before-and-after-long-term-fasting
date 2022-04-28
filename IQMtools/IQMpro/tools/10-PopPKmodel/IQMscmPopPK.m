function [] = IQMscmPopPK(nameSpace, modeltest, modelingDatasetFile, dataheaderNLME, optionsNLME, optionsSCM)
% IQMscm: Stepwise covariate search, using forward inclusion / backward
% elimination. Decision fotr inclusion or exclusion based on objective
% function value. This function takes input arguments that are used by the
% PopPK modeling workflow and thus allows to run analytic versions of
% linear PK models. 
%
% Although the algorithm is in principle applicable to categorical
% covariates with more than 2 categories, it does not really make sense and
% might introduce to many parameters that are not significant. A better
% approach is to get the modeler to think a little more and reduce the
% number of categories to 2. Covariates with multiple categories can be
% assessed after the SCM has been done by adding them one at a time.
%
% This function here will only allow for categorical covariates with 2
% categories - not more.
%
% Centering for continuous covariates is fixed at the median of the
% covariates for SCM. Structural covariates can be centered as desired.
%
% [SYNTAX]
% [] = IQMscmPopPK(nameSpace, modeltest, modelingDatasetFile, dataheaderNLME, optionsNLME, optionsSCM)
%
% [INPUT]
% Input arguments nameSpace, modeltest, modelingDatasetFile, dataheaderNLME,
% optionsNLME are identical to the ones used and documented in
% IQMbuildPopPKModelSpace and are not repeated here. 
% 
% However, instead of a model space, only a single model is allowed to be
% defined. This is checked and an error is returned if multiple models are
% defined.
%
% optionsSCM:       Structure with optional settings for the covariate search.
%   optionsSCM.outputPath:              Path where to store the log file (default: projectPathSCM)
%   optionsSCM.N_PROCESSORS_PAR:        Number of parallel model runs (default: as specified in SETUP_PATHS_TOOLS_IQMPRO)
%   optionsSCM.N_PROCESSORS_SINGLE:     Number of processors to parallelize single run (if NONMEM and MONOLIX allow for it) (default: 1)
%   optionsSCM.p_forward:               p-value for forward inclusion step (default: 0.05)
%   optionsSCM.p_backward:              p-value for forward inclusion step (default: 0.001)
%   optionsSCM.covariateTests:          Cell-array with parameter/covariate combinations to test
%                                       First element in each sub-array is the parameter name,
%                                       following are the covariates to test.
%
%                                       Example:
%                                          {
%                                             {'EMAX','WT0','SEX','HT0','ETHN'}
%                                             {'EC50','WT0','SEX','HT0','ETHN'}
%                                             {'fS'  ,'WT0','SEX','HT0','ETHN'}
%                                           }
%
%                                       If not specified or empty, then all covariates will be
%                                       tested on all parameters.
%
% If "covariateModels" is defined in modeltest, then these covariates are
% included by default and are not subject to the SCM algorithm.
%
% [OUTPUT]
% The output is in the form of a logfile that is stored in the selected
% output folder.

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Same definitions as in IQMbuilPopPKModelSpace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define parameternames in ODE and ANALYTIC models. NONLINEAR parameters need to be at the end.
% Its the required order ... checks for ODE models will be done. MONOLIX works with an analytic template. 
% For NONMEM no analytic template can be used, since numerics horribly bad. Individual ADVANs will be
% generated automatically instead.
TemplateModels = [];
TemplateModels.ParameterNames.ODE       = {'CL', 'Vc', 'Q1', 'Vp1', 'Q2', 'Vp2', 'Fiv', 'Fabs1', 'ka', 'Tlaginput1', 'Fabs0', 'Tk0input3', 'Tlaginput3', 'Frel0', 'VMAX', 'KM'};
TemplateModels.IIVdistribution.ODE      = { 'L',  'L',  'L',   'L',  'L',   'L',   'L',     'L',  'L',          'L',     'L',         'L',          'L',     'L',    'L',  'L'};
TemplateModels.ParameterNames.ANALYTIC  = {'CL', 'Vc', 'Q1', 'Vp1', 'Q2', 'Vp2', 'Fiv', 'Fabs1', 'ka', 'Tlaginput1', 'Fabs0', 'Tk0input3', 'Tlaginput3', 'Frel0'};
TemplateModels.IIVdistribution.ANALYTIC = { 'L',  'L',  'L',   'L',  'L',   'L',   'L',     'L',  'L',          'L',     'L',         'L',          'L',     'L'};
TemplateModels.Model.MONOLIX.ANALYTIC            = 'template_popPK_model_ANALYTIC_MLXTRAN.txt';
TemplateModels.Model.MONOLIX.ANALYTIC_SEQ01ABS   = 'template_popPK_model_ANALYTIC_SEQ01ABS_MLXTRAN.txt';
TemplateModels.Model.ODE                = 'template_popPK_model.txt';
TemplateModels.Model.DOSING             = 'template_popPK_dosing.dos';
TemplateModels.Model.DOSING_SEQ01ABS    = 'template_popPK_dosing_SEQ01ABS.dos';

% TemplateModels = [];
% TemplateModels.ParameterNames.ODE       = {'CL', 'Vc', 'Q1', 'Vp1', 'Q2', 'Vp2', 'Fiv', 'Fabs1', 'ka', 'Tlaginput1', 'Fabs0', 'Tk0input3', 'Tlaginput3', 'Frel0', 'VMAX', 'KM'};
% TemplateModels.IIVdistribution.ODE      = { 'L',  'L',  'L',   'L',  'L',   'L',   'L',     'L',  'L',          'L',     'L',         'L',          'L',     'L',    'L',  'L'};
% TemplateModels.ParameterNames.ANALYTIC  = {'CL', 'Vc', 'Q1', 'Vp1', 'Q2', 'Vp2', 'Fiv', 'Fabs1', 'ka', 'Tlaginput1', 'Fabs0', 'Tk0input3', 'Tlaginput3', 'Frel0'};
% TemplateModels.IIVdistribution.ANALYTIC = { 'L',  'L',  'L',   'L',  'L',   'L',   'L',     'L',  'L',          'L',     'L',         'L',          'L',     'L'};
% TemplateModels.Model.MONOLIX.ANALYTIC   = 'template_popPK_model_ANALYTIC_MLXTRAN.txt';
% TemplateModels.Model.ODE                = 'template_popPK_model.txt';
% TemplateModels.Model.DOSING             = 'template_popPK_dosing.dos';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First step: need to get all information to basically reuse IQMscm functionality
% TOOL,projectPathSCM,optionsSCM,model,dosing,data,options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get modeltest elements
try numberCompartments          = modeltest.numberCompartments;         catch, error('modeltest.numberCompartments undefined');                             end
try errorModels                 = modeltest.errorModels;                catch, error('modeltest.errorModels undefined');                                    end
try saturableClearance          = modeltest.saturableClearance;         catch, error('modeltest.saturableClearance undefined');                             end
try FACTOR_UNITS                = modeltest.FACTOR_UNITS;               catch, error('modeltest.FACTOR_UNITS undefined');                                   end
try POPvalues0                  = modeltest.POPvalues0;                 catch, error('modeltest.POPvalues0 undefined');                                     end
try POPestimate                 = modeltest.POPestimate;                catch, error('modeltest.POPestimate undefined');                                    end
try IIVestimate                 = modeltest.IIVestimate;                catch, error('modeltest.IIVestimate undefined');                                    end

try absorptionModel             = modeltest.absorptionModel;            catch, absorptionModel = 1;                                                         end % default: 1st order absorption
try lagTime                     = modeltest.lagTime;                    catch, lagTime = 0;                                                                 end % default: no lag time

try IIVvalues0                  = modeltest.IIVvalues0;                 catch, IIVvalues0           = ones(1,length(TemplateModels.ParameterNames.ODE));    end
try covarianceModel             = modeltest.covarianceModels;           catch, covarianceModel      = '';                          	                        end
try IIVdistribution             = modeltest.IIVdistribution;            catch, IIVdistribution      = {};                          	                        end
try errorParam0                 = modeltest.errorParam0;                catch, errorParam0          = [];                         	                        end
try SILENT                      = optionsNLME.SILENT;                   catch, SILENT               = 1;                         	                        end
try algorithm                   = optionsNLME.algorithm;                catch, algorithm            = [];                         	                        end
try covariateModel              = modeltest.covariateModels;            catch, covariateModel       = '';                        	                        end
try covariateModelValues        = modeltest.covariateModelValues;       catch, covariateModelValues = {};                                                   end
try COVestimate                 = modeltest.COVestimate;                catch, COVestimate          = {};                                                   end
try TOOL                        = optionsNLME.parameterEstimationTool; 	catch, TOOL                 = 'MONOLIX';                                            end

try covariateModelTV            = modeltest.covariateModelsTV;          catch, covariateModelTV     = '';                                                   end

% Check if model space is single model
multiple = 0;
if length(numberCompartments) > 1, multiple = 1; end
if iscell(errorModels), if length(errorModels) > 1, multiple = 1; else errorModels = errorModels{1}; end; end
if iscell(errorParam0), if length(errorParam0) > 1, multiple = 1; else errorParam0 = errorParam0{1}; end; end
if length(saturableClearance) > 1, multiple = 1; end
if length(absorptionModel) > 1, multiple = 1; end
if length(lagTime) > 1, multiple = 1; end
if iscell(POPvalues0), if length(POPvalues0)>1, multiple = 1; else POPvalues0 = POPvalues0{1}; end; end
if iscell(POPestimate), if length(POPestimate)>1, multiple = 1; else POPestimate = POPestimate{1}; end; end
if iscell(IIVestimate), if length(IIVestimate)>1, multiple = 1; else IIVestimate = IIVestimate{1}; end; end
if iscell(IIVvalues0), if length(IIVvalues0)>1, multiple = 1; else IIVvalues0 = IIVvalues0{1}; end; end
if iscell(covarianceModel), if length(covarianceModel)>1, multiple = 1; else covarianceModel = covarianceModel{1}; end; end
if iscell(covariateModel), if length(covariateModel)>1, multiple = 1; else covariateModel = covariateModel{1}; end; end
if multiple,
    error('Modelspace is larger than a single model.');
end

% Get TOOL
if strcmpi(TOOL,'nonmem'),
    FLAG_NONMEM                 = 1;
elseif strcmpi(TOOL,'monolix'),
    FLAG_NONMEM                 = 0;
else
    error('Unknown definition of tool in optionsNLME.parameterEstimationTool.');
end

% Get projectPathSCM
projectPathSCM              	= ['../Models/' nameSpace];

% Get optionsSCM
%   Already defined as input argument

% Get model
%   Not needed ... will be handled differently

% Get dosing
%   Not needed ... will be handled differently

% Get data
[xdummyx,f,e]                   = fileparts(modelingDatasetFile);
data                            = [];
data.dataFileName                = [f e];
data.header                     = dataheaderNLME;
data.dataRelPathFromProject     = '../../Data';

% Get options for single model generation by copying
options                         = [];
options.POPestimate             = POPestimate;
options.POPvalues0              = POPvalues0;
options.IIVdistribution         = IIVdistribution;
options.IIVestimate             = IIVestimate;
options.IIVvalues0              = IIVvalues0;
options.errorModels             = errorModels;
options.errorParam0             = errorParam0;
options.covarianceModel         = covarianceModel;
options.SILENT                  = SILENT;
options.algorithm               = algorithm;
options.covariateModel          = covariateModel;
options.covariateModelValues    = covariateModelValues;
options.COVestimate             = COVestimate;

options.covariateModelTV        = covariateModelTV;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Second step: Do this IQMscm but with different model generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Handle optionsSCM
try p_forward               = optionsSCM.p_forward;                 catch, p_forward                = 0.05;                     end
try p_backward              = optionsSCM.p_backward;                catch, p_backward               = 0.001;                    end
try covariateTests          = optionsSCM.covariateTests;            catch, covariateTests           = {};                       end
try outputPath              = optionsSCM.outputPath;                catch, outputPath               = projectPathSCM;           end
try N_PROCESSORS_PAR        = optionsSCM.N_PROCESSORS_PAR;          catch, N_PROCESSORS_PAR         = getN_PROCESSORS_PARIQM(); end
try N_PROCESSORS_SINGLE     = optionsSCM.N_PROCESSORS_SINGLE;       catch, N_PROCESSORS_SINGLE      = 1;                        end

% Get previously defined covariate model (fixed covariate model)
try covariateModel          = options.covariateModel;               catch, covariateModel           = '';                       end
try covariateModelValues    = options.covariateModelValues;         catch, covariateModelValues     = {};                       end
try COVestimate             = options.COVestimate;                  catch, COVestimate              = {};                       end

% Check if NONMEM and SAEM ... then show error
if strcmpi(TOOL,'NONMEM') && strcmpi(options.algorithm.METHOD,'SAEM'),
    error('Stepwise covariate search with NONMEM SAEM is not a good idea. Please use FO(CE(I)) or MONOLIX.');
end

% Create the projectFolder
try rmdir(projectPathSCM,'s'); catch end
warning off
mkdir(projectPathSCM);
warning on

% Load dataset to assess number of categorical elements
oldpath = pwd();
cd(projectPathSCM);
dataContents = IQMloadCSVdataset([data.dataRelPathFromProject '/' data.dataFileName]);
cd(oldpath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check unsupported administration routes
% ADM=1 is covering 1st order, 0th order, and mixed absorption for the popPK
% workflow. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
ADMunique = unique(dataContents.ADM);
if sum(ADMunique>2)~=0,
    error('Unknown/unsupported ADM entries in dataset. 0: observation, 1: 0th/1st/mixed order absorption, 2: IV (bolus or infusion).');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check lag time and absorption model settings.
% If no absorption doses (ADM=1) present in the dataset and lagTime and 
% absorptionModel variables not on default values, then reset them to 
% default values. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
% Check if absorption dose data present (ADM=1)
FLAG_ABSORPTION_DATA_PRESENT = 1;
if ~ismember(1,ADMunique),
    FLAG_ABSORPTION_DATA_PRESENT = 0;
    if sum(lagTime) ~= 0,
        warning('Dataset does not contain absorption doses (ADM=1). Estimation of lag time not considered.');
        lagTime = 0;
    end
    if sum(absorptionModel) ~= 0,
        warning('Dataset does not contain absorption doses (ADM=1). No absorption model considered.');
        absorptionModel = 1;
    end
end

%% Intermediate step: Check NONMEM issue 
% Fiv not possible to estimate with NONMEM since RATE not handled together
% with bioavailability. In this case MONOLIX needs to be used.
% Only needs to be checked if IV data present in the dataset.

FLAG_NONMEM__RATE_BIOAVAILABILITY_ISSUE__OK = 1;

if ~isempty(find(dataContents.ADM==2)),
    % IV data present ... need to check
    
    % Get index Fiv
    ix_Fiv = strmatchIQM('Fiv',TemplateModels.ParameterNames.ODE,'exact');
    
    % If Fiv initial guess not 1 then NONMEM not ok
    if POPvalues0(ix_Fiv) ~= 1,
        FLAG_NONMEM__RATE_BIOAVAILABILITY_ISSUE__OK = 0;
    end
    
    % If Fiv estimated then NONMEM not ok
    if POPestimate(ix_Fiv) ~= 0,
        FLAG_NONMEM__RATE_BIOAVAILABILITY_ISSUE__OK = 0;
    end
    
    % If FACTOR_UNITS not equal 1 then NONMEM not ok
    if FACTOR_UNITS~=1,
        FLAG_NONMEM__RATE_BIOAVAILABILITY_ISSUE__OK = 0;
    end
    
    % Error message if needed
    if FLAG_NONMEM && ~FLAG_NONMEM__RATE_BIOAVAILABILITY_ISSUE__OK,
        error(sprintf('NONMEM can not be used with these settings (estimation of Fiv or non-unity initial guess for Fiv or FACTOR_UNITS not equal to 1).\nNONMEM is unable to combine bioavailability parameter with infusion RATE. Simple fix: For these model\nsettings, please use the more modern software MONOLIX, which can do this correctly without a lot of work-arounds.'));
    end
end

%% Change data path to match the nested model folder structure
data.dataRelPathFromProject = ['../' data.dataRelPathFromProject];

% If NONMEM used then check if IMP set to 1
if strcmpi(TOOL,'nonmem') && strcmpi(options.algorithm.METHOD,'saem') && options.algorithm.IMPORTANCESAMPLING == 0,
    error('When using NONMEM/SAEM, please set the options.algorithm.IMPORTANCESAMPLING=1.');
end

% Create the BASE model (MODEL_BASE) and load header
BASEmodelFolder             = [projectPathSCM '/MODEL_BASE'];
dataRelPathFromProjectPath  = '../../../Data';
createPopPK_NLMEproject_ODE_Analytic_IQM(FLAG_NONMEM,FLAG_ABSORPTION_DATA_PRESENT,'MODEL_BASE',TemplateModels,FACTOR_UNITS, ...
            numberCompartments,saturableClearance, ...
            absorptionModel,lagTime, ...
            data,BASEmodelFolder,dataRelPathFromProjectPath,options, ...
            FLAG_NONMEM__RATE_BIOAVAILABILITY_ISSUE__OK, dataContents);
projectInfo                 = parseNLMEprojectHeaderIQM(BASEmodelFolder);
                                       
% Get parameter names, continuous and categorical covariate names from the
% models header
parNames            = projectInfo.PARAMNAMES;
covNames            = projectInfo.COVNAMES;
catNames            = projectInfo.CATNAMES;
covcatNames         = [covNames catNames];

% If list of covariates empty then generate it - to include all covariates
% on all parameters
if isempty(covariateTests)
    covariateTests = cell(length(parNames),1);
    for k1=1:length(parNames),
        cTk{1} = parNames{k1};
        for k2=1:length(covcatNames),
            cTk{k2+1} = covcatNames{k2};
        end
        covariateTests{k1} = cTk(:)';
    end
end

% Expand the covariat settings to matching lists for parameters and covariate names
covSearch                   = [];
covSearch.paramNamesTest    = {};
covSearch.covcatNamesTest   = {};
for k1=1:length(covariateTests),
    for k2=2:length(covariateTests{k1}),
        covSearch.paramNamesTest{end+1}  = covariateTests{k1}{1};
        covSearch.covcatNamesTest{end+1} = covariateTests{k1}{k2};
    end
end

% Check lists if selected covariates and parameters present in model and dataset
for k=1:length(covSearch.paramNamesTest),
    ix = strmatchIQM(covSearch.paramNamesTest{k},parNames,'exact');
    if isempty(ix),
        error('Parameter %s is not estimated in the model.',covSearch.paramNamesTest{k});
    end
end

for k=1:length(covSearch.covcatNamesTest),
    ix = strmatchIQM(covSearch.covcatNamesTest{k},covcatNames,'exact');
    if isempty(ix),
        error('Covariate %s is not present in the dataset.',covSearch.covcatNamesTest{k});
    end
end

% Determine structural covariate information (structural covariates are
% covariates that already are in the model by options.covariateModel)
covStructural               = [];
covStructural.paramNames    = {};
covStructural.covcatNames   = {};
covStructural.value         = [];
covStructural.estimate      = [];
terms = strrep(strrep(explodePCIQM(covariateModel,',','{','}'),'{',''),'}','');
for k=1:length(terms),
    x = explodePCIQM(terms{k});
    for k2=2:length(x),
        covStructural.paramNames{end+1} = x{1};
        covStructural.covcatNames{end+1} = x{k2};
        if isempty(covariateModelValues),
            covStructural.value(end+1) = 0.1;
        else
            covStructural.value(end+1) = covariateModelValues{k}(k2-1);
        end
        if isempty(COVestimate),
            covStructural.estimate(end+1) = 1;
        else
            covStructural.estimate(end+1) = COVestimate{k}(k2-1);
        end
    end
    
end

% Remove structural covariates from covSearch - since these are already in
% the model
ix_remove = [];
for k=1:length(covStructural.paramNames),
    par = covStructural.paramNames{k};
    cov = covStructural.covcatNames{k};
    ix1 = strmatchIQM(par,covSearch.paramNamesTest,'exact');
    ix2 = strmatchIQM(cov,covSearch.covcatNamesTest,'exact');
    if ~isempty(intersect(ix1,ix2)),
        ix_remove(end+1) = intersect(ix1,ix2);
    end
end
if ~isempty(ix_remove),
    covSearch.paramNamesTest(ix_remove) = [];
    covSearch.covcatNamesTest(ix_remove) = [];
end

% Determine delta Objective functions required for covariates
% DF: 1 for continuous, N-1 for categorical - since N=2 for categorical ...
% it is the same.
for k=1:length(covSearch.covcatNamesTest),
    % Check if categorical
    if ~isempty(strmatchIQM(covSearch.covcatNamesTest{k},covNames,'exact')),
        % Continuous covariate
        N_DF = 1;
    else
        % Categorical covariate
        N_categories = length(unique(dataContents.(covSearch.covcatNamesTest{k})));
        if N_categories > 2,
            error('Categorical covariate "%s" has more than 2 categories - this is not allowed.',covSearch.covcatNamesTest{k});
        end
        if N_categories < 2,
            error('Categorical covariate "%s" has less than 2 categories - this is not allowed.',covSearch.covcatNamesTest{k});
        end
        N_DF = 1;        
    end
    covSearch.delta_forward(k)   = chi2invIQM(1-p_forward,N_DF);
    covSearch.delta_backward(k)  = chi2invIQM(1-p_backward,N_DF);
end

% Save information for later
covSearchBackUp = covSearch;

% Open report file
warning off
mkdir(outputPath);
warning on
fid = fopen([outputPath '/SCMlogfile.txt'],'w');

fprintf(fid,'****************************************************************************\n');
fprintf(fid,'* Covariate search using forward inclusion and backward elimination method *\n');
fprintf(fid,'****************************************************************************\n');
fprintf(fid,'Username: %s\n',usernameIQM());
fprintf(fid,'Date:     %s\n',datestr(now,'yyyy-mmm-DD HH:MM'));
fprintf(fid,'****************************************************************************\n');
fprintf(fid,'\n');
fprintf(fid,'p-value forward:      %g\n',p_forward);
fprintf(fid,'p-value backward:     %g\n',p_backward);
fprintf(fid,'\n');
fprintf(fid,'Decisions based on delta objective function.\n');
fprintf(fid,'\n');

fprintf(fid,'************************************\n');
fprintf(fid,'* RUNNING BASE MODEL (MODEL_BASE)  *\n');
fprintf(fid,'************************************\n');
fprintf(fid,'\n');

% Run the BASE model (MODEL_BASE)
BASEmodelFolder = [projectPathSCM '/MODEL_BASE'];
IQMrunNLMEproject(BASEmodelFolder,N_PROCESSORS_SINGLE);
RESULTS         = parseNLMEprojectResults(BASEmodelFolder);

% Get OFV and optimal parameters for BASE model
BASE_OFV            = RESULTS.objectivefunction.OBJ;

% If NaN then fit crashed ... pbly NONMEM error
if isnan(BASE_OFV),
    error('Base model crashed ... ');
end

% Report OFV for BASE model
fprintf(fid,'Objective function value for BASE model: %g\n',BASE_OFV);
if strcmpi(TOOL,'NONMEM'),
    fprintf(fid,'\n');
    fprintf(fid,'%s\n',RESULTS.termination_info{1});
end
fprintf(fid,'\n');

% Update model options with optimal BASE model parameters
options.POPvalues0  = RESULTS.rawParameterInfo.fixedEffects.values;
options.IIVvalues0  = RESULTS.rawParameterInfo.randomEffects.values;
options.errorParam0 = RESULTS.rawParameterInfo.errorParameter.values;

% If sequential 0/1 order absorption need to add Tlaginput 1 in between 
if absorptionModel == 3,
    POPvalues0 = NaN(1,length(options.POPvalues0)+1);
    POPvalues0([1:9 11:end]) = options.POPvalues0;
    options.POPvalues0 = POPvalues0;
    
    IIVvalues0 = NaN(1,length(options.IIVvalues0)+1);
    IIVvalues0([1:9 11:end]) = options.IIVvalues0;
    options.IIVvalues0 = IIVvalues0;
end

% Update data path again for the level models
data.dataRelPathFromProject = ['../' data.dataRelPathFromProject];

% Forward inclusion
fprintf(fid,'************************************\n');
fprintf(fid,'* FORWARD INCLUSION                *\n');
fprintf(fid,'************************************\n');
fprintf(fid,'\n');

OLD_OFV                     = BASE_OFV;
covariatesForward           = {};
covariateModelValuesForward = {};
COVestimateForward          = {};
countModel                  = 1;
Nmaxmodels                  = (length(covSearch.covcatNamesTest)^2-length(covSearch.covcatNamesTest))/2+length(covSearch.covcatNamesTest);
continueForwardSearch       = 1;
levelCount                  = 1;
levelModel                  = '';
ForwardModelName            = BASEmodelFolder;
ModelNameAll                = {};

while continueForwardSearch,
    
    % Run through remaining parameter/covariate combinations
    % Generate all models first and then run all at once
    covariatesTestForward           = {};
    covariateModelValuesTESTForward = {};
    COVestimateTESTForward          = {};
    covariateModelAll               = {};
    
    for k1=1:length(covSearch.paramNamesTest),
        
        % Reset covariate setting to previous level
        covariatesTest          = covariatesForward;
        covariateModelValues    = covariateModelValuesForward;
        COVestimate             = COVestimateForward;
        
        % Get new parameter/covariate combination to test
        param = covSearch.paramNamesTest{k1};
        cov   = covSearch.covcatNamesTest{k1};
        delta = covSearch.delta_forward(k1);
        value = 0.1; % Setting default value to 0.1 (NONMEM does not allow 0)
        
        % Build covariate test structure and covariateModelValues and COVestimate
        ix = [];
        for k2=1:length(covariatesTest),
            if strcmp(covariatesTest{k2}{1},param),
                ix = k2;
            end
        end
        if isempty(ix),
            covariatesTest{end+1}{1}        = param;
            covariatesTest{end}{2}          = cov;
            covariateModelValues{end+1}     = value;
            COVestimate{end+1}              = 1;
        else
            covariatesTest{ix}{end+1}       = cov;
            covariateModelValues{ix}(end+1) = value;
            COVestimate{ix}(end+1)          = 1;
        end
        covariatesTest                      = covariatesTest(:);
        
        % Save for later
        covariatesTestForward{k1}           = covariatesTest;
        covariateModelValuesTESTForward{k1} = covariateModelValues;
        COVestimateTESTForward{k1}          = COVestimate;
        
        % Convert to required format
        covariateModel                      = '';
        for k3=1:size(covariatesTest,1),
            if ~isempty(covariatesTest{k3}),
                text                        = sprintf('%s,',covariatesTest{k3}{:});
                covariateModel              = [covariateModel '{' text(1:end-1) '},'];
            end
        end
        covariateModel                      = covariateModel(1:end-1);
        covariateModelAll{k1}               = covariateModel;
        
        % Determine COVestimate and covariateModelValues
        COVestimate = {};
        covariateModelValues = {};
        for k3=1:size(covariatesTest,1),
            COVestimate{k3}                 = ones(1,length(covariatesTest{k3})-1);
            covariateModelValues{k3}        = 0.1*ones(1,length(covariatesTest{k3})-1);
        end
        
        % Add structural covariates to covariate model
        for k3=1:length(covStructural.paramNames),
            param    = covStructural.paramNames{k3};
            cov      = covStructural.covcatNames{k3};
            value    = covStructural.value(k3);
            estimate = covStructural.estimate(k3);
            % Find where to add
            ix = [];
            for k4=1:size(covariatesTest,1),
                if strcmp(param,covariatesTest{k4}{1}),
                    ix = k4;
                    break;
                end
            end
            if ~isempty(ix),
                % add info
                covariatesTest{ix}{end+1} = cov;
                covariateModelValues{ix}(end+1) = value;
                COVestimate{ix}(end+1) = estimate;
            else
                % new param
                covariatesTest{end+1}{1} = param;
                covariatesTest{end}{2} = cov;
                covariateModelValues{end+1}(1) = value;
                COVestimate{end+1}(1) = estimate;
            end
            covariatesTest = covariatesTest(:);
        end
        
        % Convert to required format for the search things
        covariateModel                      = '';
        for k3=1:size(covariatesTest,1),
            if ~isempty(covariatesTest{k3}),
                text                        = sprintf('%s,',covariatesTest{k3}{:});
                covariateModel              = [covariateModel '{' text(1:end-1) '},'];
            end
        end
        covariateModel                      = covariateModel(1:end-1);
        
        % Set covariateModel
        options.covariateModel              = covariateModel;
        options.covariateModelValues        = covariateModelValues;
        options.COVestimate                 = COVestimate;
        
        % Create model
        ModelName                           = sprintf('MODEL_%s',preFillCharIQM(countModel,length(num2str(Nmaxmodels)),'0'));
        FolderName                          = [projectPathSCM sprintf('/FW_LEVEL_%d',levelCount)];
        ModelFolder                         = [FolderName '/' ModelName];
        ModelNameAll{k1}                    = ModelFolder;
        
        % For popPK all parameters need to be positive. 0 is not allowed due to L or maybe G transformation
        % Here we fix 0 values to 1e-10 if present
        options.POPvalues0(options.POPvalues0==0) = 1e-10;
        
        % Create project    
        dataRelPathFromProjectPath = '../../../../Data';
        createPopPK_NLMEproject_ODE_Analytic_IQM(FLAG_NONMEM,FLAG_ABSORPTION_DATA_PRESENT,ModelName,TemplateModels,FACTOR_UNITS, ...
            numberCompartments,saturableClearance,...
            absorptionModel,lagTime,...
            data,ModelFolder,dataRelPathFromProjectPath,options,...
            FLAG_NONMEM__RATE_BIOAVAILABILITY_ISSUE__OK, dataContents);
        
        countModel                          = countModel + 1;
    end
    
    % Run the level models - do not create GoF plots
    IQMrunNLMEprojectFolder(FolderName,N_PROCESSORS_PAR,N_PROCESSORS_SINGLE,1);
    
    % Read the results 
    NLME_ORDER_CRITERION    = ''; % Do not order them
    LEVEL_RESULTS           = parseProjectFolderResultsIQM(FolderName,NLME_ORDER_CRITERION);
    IQMfitsummaryAll(FolderName,FolderName,NLME_ORDER_CRITERION);
    
    % Determine significance 
    LEVEL_SIGNIFICANCE      = [];
    % Get OBJ
    LEVEL_OBJ               = [LEVEL_RESULTS.OBJ];
    
    % Create a table with output information    
    if strcmpi(TOOL,'NONMEM'),
        tableText               = cell(1,10);
        tableText(1,1:2)        = {'<TT>' sprintf('Level %d: %s',levelCount,levelModel)};
        tableText(2,:)          = {'<TH>','Model','Covariate tested','OBJ(cov)','OBJ(prev)','OBJ(prev)-OBJ(cov)','Significance','MINIMIZATION SUCCESSFUL','MINIMIZATION TERMINATED','PROBLEMS'};
    else
        tableText               = cell(1,7);
        tableText(1,1:2)        = {'<TT>' sprintf('Level %d: %s',levelCount,levelModel)};
        tableText(2,:)          = {'<TH>','Model','Covariate tested','OBJ(cov)','OBJ(prev)','OBJ(prev)-OBJ(cov)','Significance'};
    end
    for k1=1:length(covSearch.paramNamesTest),
        tableText{k1+2,1}   = '<TR>';
        tableText{k1+2,2}   = LEVEL_RESULTS(k1).model;
        tableText{k1+2,3}   = sprintf('%s on %s',covSearch.covcatNamesTest{k1},covSearch.paramNamesTest{k1});
        tableText{k1+2,4}   = round(100*LEVEL_OBJ(k1))/100;
        tableText{k1+2,5}   = round(100*OLD_OFV)/100;
        tableText{k1+2,6}   = round(100*(OLD_OFV-LEVEL_OBJ(k1)))/100;
        if OLD_OFV-LEVEL_OBJ(k1) >= covSearch.delta_forward(k1),
            LEVEL_SIGNIFICANCE(k1) = 1;
            tableText{k1+2,7}   = 'YES';
        else
            LEVEL_SIGNIFICANCE(k1) = 0;
            tableText{k1+2,7}   = '-';
        end
        if strcmpi(TOOL,'NONMEM'),
            if ~isempty(strfind(LEVEL_RESULTS(k1).termination_info,'MINIMIZATION SUCCESSFUL')),
                tableText{k1+2,8} = 'YES';
            else
                tableText{k1+2,8} = '-';
            end            
            if ~isempty(strfind(LEVEL_RESULTS(k1).termination_info,'MINIMIZATION TERMINATED')),
                tableText{k1+2,9} = 'YES';
            else
                tableText{k1+2,9} = '-';
            end     
            if ~isempty(strfind(LEVEL_RESULTS(k1).termination_info,'PROBLEMS')),
                tableText{k1+2,10} = 'YES';
            else
                tableText{k1+2,10} = '-';
            end            
        end
    end        
    
    % Determine which covariate leads to the largest drop in OBJ
    [~,ixBEST] = max(OLD_OFV-LEVEL_OBJ);
    
    % If none significant then stop forward inclusion here
    if sum(LEVEL_SIGNIFICANCE) == 0,
        % Nothing significant anymore - stop forward search
        continueForwardSearch = 0;
        
        % Finalize table for this level by adding footer
        tableText(end+1,1:2) = {'<TF>' 'No significant decrease in objective function => End of forward inclusion.'};
    else
        % Finalize table for this level by adding footer
        tableText(end+1,1:2) = {'<TF>' sprintf('Retained covariate model for next level: %s (%s on %s)\n',LEVEL_RESULTS(ixBEST).model,covSearch.covcatNamesTest{ixBEST},covSearch.paramNamesTest{ixBEST})};

        % Setup for new level
        % -------------------
        
        OLD_OFV                             = LEVEL_OBJ(ixBEST);
        
        % Update options with optimal model parameters from last run
        options.POPvalues0                  = LEVEL_RESULTS(ixBEST).rawParameterInfo.fixedEffects.values;
        options.IIVvalues0                  = LEVEL_RESULTS(ixBEST).rawParameterInfo.randomEffects.values;
        options.errorParam0                 = LEVEL_RESULTS(ixBEST).rawParameterInfo.errorParameter.values;

        % If sequential 0/1 order absorption need to add Tlaginput 1 in between
        if absorptionModel == 3,
            POPvalues0 = NaN(1,length(options.POPvalues0)+1);
            POPvalues0([1:9 11:end]) = options.POPvalues0;
            options.POPvalues0 = POPvalues0;
            
            IIVvalues0 = NaN(1,length(options.IIVvalues0)+1);
            IIVvalues0([1:9 11:end]) = options.IIVvalues0;
            options.IIVvalues0 = IIVvalues0;
        end
           
        % Remove covariate from search
        covSearch.covcatNamesTest(ixBEST)   = [];
        covSearch.delta_forward(ixBEST)     = [];
        covSearch.delta_backward(ixBEST)    = [];
        covSearch.paramNamesTest(ixBEST)    = [];
        
        % Set new base covariates
        covariatesForward                   = covariatesTestForward{ixBEST};
        covariateModelValuesForward         = covariateModelValuesTESTForward{ixBEST};
        COVestimateForward                  = COVestimateTESTForward{ixBEST};
        
        % Update estimated covariate coefficient for next level - with all new estimates!
        covInfo = LEVEL_RESULTS(ixBEST).rawParameterInfo.covariate;
        
        for k1x = 1:length(covariatesForward),
            paramUpdate = covariatesForward{k1x}{1};
            for k2x = 2:length(covariatesForward{k1x}),
                covUpdate = covariatesForward{k1x}{k2x};
                % Find parameter and covariate
                matchIX = [];
                for kkkx=1:length(covInfo.names),
                    ixParam = strfind(covInfo.names{kkkx},paramUpdate);
                    ixCov = strfind(covInfo.names{kkkx},covUpdate);
                    if ~isempty(ixParam) && ~isempty(ixCov),
                        matchIX(end+1) = kkkx;
                    end
                end
                if length(matchIX) ~= 1,
                    error('Problem with getting covariate coefficient value.');
                end
                % Get the estimated value
                value = covInfo.values(matchIX);
                % Add value to covariateModelValuesForward
                ix = [];
                for k2=1:length(covariatesForward),
                    if strcmp(covariatesForward{k2}{1},paramUpdate),
                        ix = k2;
                    end
                end
                covariateModelValuesForward{ix}(strmatchIQM(covUpdate,covariatesForward{ix},'exact')-1) = value;
            end
        end
        
        levelModel          = covariateModelAll{ixBEST};
        ForwardModelName    = ModelNameAll{ixBEST};
        levelCount          = levelCount+1;
    end
    
    % Export table for the level into the log file
    textDisplay = IQMconvertCellTable2ReportTable(tableText,'text');
    fprintf(fid,'%s\n',textDisplay);
end

% Get forward results
FORWARD_OFV = OLD_OFV;
forwardModel = levelModel;

% Report forward model
fprintf(fid,'************************************\n');
fprintf(fid,'* FORWARD MODEL RESULTS            *\n');
fprintf(fid,'************************************\n');
fprintf(fid,'\n');

fprintf(fid,'Forward model:             %s\n',ForwardModelName);
fprintf(fid,'Covariates:                %s\n',forwardModel);
fprintf(fid,'Objective function value:  %g\n',FORWARD_OFV);
fprintf(fid,'\n');


% Backward elimination
fprintf(fid,'************************************\n');
fprintf(fid,'* BACKWARD ELIMINATION             *\n');
fprintf(fid,'************************************\n');
fprintf(fid,'\n');

% Create covSearch structure for backward search
covSearchBackward = covSearchBackUp;
ix_keep = [];
for k1=1:size(covariatesForward),
    param = covariatesForward{k1}{1};
    for k2=2:length(covariatesForward{k1}),
        cov = covariatesForward{k1}{k2};
        ix_keep(end+1) = intersect(strmatchIQM(param,covSearchBackward.paramNamesTest,'exact'),strmatchIQM(cov,covSearchBackward.covcatNamesTest,'exact'));
    end
end
covSearchBackward.covcatNamesTest = covSearchBackward.covcatNamesTest(ix_keep);
covSearchBackward.paramNamesTest = covSearchBackward.paramNamesTest(ix_keep);
covSearchBackward.delta_forward = covSearchBackward.delta_forward(ix_keep);
covSearchBackward.delta_backward = covSearchBackward.delta_backward(ix_keep);

% Add parameter values
for k1=1:length(covSearchBackward.paramNamesTest),
    paramName = covSearchBackward.paramNamesTest{k1};
    covName   = covSearchBackward.covcatNamesTest{k1};
    for k2=1:size(covariatesForward),
        if strcmp(paramName,covariatesForward{k2}{1}),
            ix = k2;
            break;
        end
    end
    ix2 = strmatchIQM(covName,covariatesForward{ix}(2:end),'exact');
    value = covariateModelValuesForward{ix}(ix2);
    covSearchBackward.value(k1) = value;
end

% Do the search
covariateModelAll  = {};
levelCount         = 1;
countModel         = 1;
BackwardModelName  = ForwardModelName;
ModelNameAll       = {};

% Check if backward search needed
if ~isempty(covariatesForward),
    continueBackwardSearch = 1;
else
    continueBackwardSearch = 0;
end

% OBJ:  Start with FORWARD_OFV and remove one by one ... until deltaOFV>... for
while continueBackwardSearch,
    
    % Run through remaining parameter/covariate combinations
    % Generate all models first and then run all at once
    
    for k1=1:length(covSearchBackward.paramNamesTest),
        
        % Reset covariate setting to previous level
        covariatesTest  = covSearchBackward;
        
        % Remove parameter to test
        param                   = covariatesTest.paramNamesTest;
        cov                     = covariatesTest.covcatNamesTest;
        value                   = covariatesTest.value;
        param(k1)               = [];
        cov(k1)                 = [];
        value(k1)               = [];
        
        % Build covariate test structure
        x                       = unique(param);
        covTestStructure        = {};
        covariateModelValues    = cell(1,length(x));
        COVestimate             = cell(1,length(x));
        
        for k2=1:length(x),
            covTestStructure{k2}{1}         = x{k2};
        end
        covTestStructure                    = covTestStructure(:);
        for k2=1:length(param),
            ix = strmatchIQM(param{k2},x,'exact');
            covTestStructure{ix}{end+1} 	= cov{k2};
            covariateModelValues{ix}(end+1) = value(k2);
            COVestimate{ix}(end+1)          = 1;
        end
        
        % Convert to required format
        covariateModel                      = '';
        for k3=1:size(covTestStructure,1),
            if ~isempty(covTestStructure{k3}),
                text                        = sprintf('%s,',covTestStructure{k3}{:});
                covariateModel              = [covariateModel '{' text(1:end-1) '},'];
            end
        end
        covariateModel                      = covariateModel(1:end-1);
        covariateModelAll{k1}               = covariateModel;
        
        % Add structural covariates to covariate model
        for k3=1:length(covStructural.paramNames),
            param                           = covStructural.paramNames{k3};
            cov                             = covStructural.covcatNames{k3};
            value                           = covStructural.value(k3);
            estimate                        = covStructural.estimate(k3);
            % Find where to add
            ix = [];
            for k4=1:size(covTestStructure,1),
                if strcmp(param,covTestStructure{k4}{1}),
                    ix = k4;
                    break;
                end
            end
            if ~isempty(ix),
                % add info
                covTestStructure{ix}{end+1}         = cov;
                covariateModelValues{ix}(end+1)     = value;
                COVestimate{ix}(end+1)              = estimate;
            else
                % new param
                covTestStructure{end+1}{1}          = param;
                covTestStructure{end}{2}            = cov;
                covariateModelValues{end+1}(1)      = value;
                COVestimate{end+1}(1)               = estimate;
            end
            covTestStructure                        = covTestStructure(:);
        end
        
        % Convert to required format for the search things
        covariateModel                      = '';
        for k3=1:size(covTestStructure,1),
            if ~isempty(covTestStructure{k3}),
                text                        = sprintf('%s,',covTestStructure{k3}{:});
                covariateModel              = [covariateModel '{' text(1:end-1) '},'];
            end
        end
        covariateModel                      = covariateModel(1:end-1);
        
        % Set covariateModel
        options.covariateModel              = covariateModel;
        options.covariateModelValues        = covariateModelValues;
        options.COVestimate                 = COVestimate;
        
        % Create model
        ModelName                           = sprintf('MODEL_%s',preFillCharIQM(countModel,length(num2str(Nmaxmodels)),'0'));
        FolderName                          = [projectPathSCM sprintf('/BW_LEVEL_%d',levelCount)];
        ModelFolder                         = [FolderName '/' ModelName];
        ModelNameAll{k1}                    = ModelFolder;
        
        % For popPK all parameters need to be positive. 0 is not allowed due to L or maybe G transformation
        % Here we fix 0 values to 1e-10 if present
        options.POPvalues0(options.POPvalues0==0) = 1e-10;
        
        % Create project    
        dataRelPathFromProjectPath = '../../../../Data';
        createPopPK_NLMEproject_ODE_Analytic_IQM(FLAG_NONMEM,FLAG_ABSORPTION_DATA_PRESENT,ModelName,TemplateModels,FACTOR_UNITS,...
            numberCompartments,saturableClearance,...
            absorptionModel,lagTime,...
            data,ModelFolder,dataRelPathFromProjectPath,options,...
            FLAG_NONMEM__RATE_BIOAVAILABILITY_ISSUE__OK, dataContents);
        
        countModel                          = countModel + 1;
    end
    
    % Run the level models - do not create GoF plots
    IQMrunNLMEprojectFolder(FolderName,N_PROCESSORS_PAR,N_PROCESSORS_SINGLE,1);
    
    % Read the results
    NLME_ORDER_CRITERION = ''; % Do not order them
    LEVEL_RESULTS = parseProjectFolderResultsIQM(FolderName,NLME_ORDER_CRITERION);
    IQMfitsummaryAll(FolderName,FolderName,NLME_ORDER_CRITERION);
    
    % Determine significance 
    LEVEL_SIGNIFICANCE = [];
    LEVEL_OBJ = [LEVEL_RESULTS.OBJ];
    
    
    % Create a table with output information    
    if strcmpi(TOOL,'NONMEM'),
        tableText               = cell(1,10);
        tableText(1,1:2)        = {'<TT>' sprintf('Level %d: %s',levelCount,levelModel)};
        tableText(2,:)          = {'<TH>','Model','Covariate removed','OBJ(cov)','OBJ(prev)','OBJ(prev)-OBJ(cov)','Significance','MINIMIZATION SUCCESSFUL','MINIMIZATION TERMINATED','PROBLEMS'};
    else
        tableText               = cell(1,7);
        tableText(1,1:2)        = {'<TT>' sprintf('Level %d: %s',levelCount,levelModel)};
        tableText(2,:)          = {'<TH>','Model','Covariate tested','OBJ(cov)','OBJ(prev)','OBJ(prev)-OBJ(cov)','Significance'};
    end
    for k1=1:length(covSearchBackward.paramNamesTest),
        tableText{k1+2,1}   = '<TR>';
        tableText{k1+2,2}   = LEVEL_RESULTS(k1).model;
        tableText{k1+2,3}   = sprintf('%s on %s',covSearchBackward.covcatNamesTest{k1},covSearchBackward.paramNamesTest{k1});
        tableText{k1+2,4}   = round(100*LEVEL_OBJ(k1))/100;
        tableText{k1+2,5}   = round(100*OLD_OFV)/100;
        tableText{k1+2,6}   = round(100*(-OLD_OFV+LEVEL_OBJ(k1)))/100;
        if -OLD_OFV+LEVEL_OBJ(k1) >= covSearchBackward.delta_backward(k1),
            LEVEL_SIGNIFICANCE(k1) = 1;
            tableText{k1+2,7}   = 'YES';
        else
            LEVEL_SIGNIFICANCE(k1) = 0;
            tableText{k1+2,7}   = '-';
        end
        if strcmpi(TOOL,'NONMEM'),
            if ~isempty(strfind(LEVEL_RESULTS(k1).termination_info,'MINIMIZATION SUCCESSFUL')),
                tableText{k1+2,8} = 'YES';
            else
                tableText{k1+2,8} = '-';
            end            
            if ~isempty(strfind(LEVEL_RESULTS(k1).termination_info,'MINIMIZATION TERMINATED')),
                tableText{k1+2,9} = 'YES';
            else
                tableText{k1+2,9} = '-';
            end     
            if ~isempty(strfind(LEVEL_RESULTS(k1).termination_info,'PROBLEMS')),
                tableText{k1+2,10} = 'YES';
            else
                tableText{k1+2,10} = '-';
            end            
        end
    end        
    
    % Determine which covariate leads to the largest increase in OBJ
    [~,ixBEST] = min(-OLD_OFV+LEVEL_OBJ);
    
    % If all significant then stop backward elimination here
    if sum(LEVEL_SIGNIFICANCE~=1) == 0,
        continueBackwardSearch = 0;
        % Finalize table for this level by adding footer
        tableText(end+1,1:2) = {'<TF>' 'No insignificant increase in objective function => End of backward elimination.'};
    else
        % Finalize table for this level by adding footer
        tableText(end+1,1:2) = {'<TF>' sprintf('Removed covariate model for next level: %s (%s on %s)\n',LEVEL_RESULTS(ixBEST).model,covSearchBackward.covcatNamesTest{ixBEST},covSearchBackward.paramNamesTest{ixBEST})};
        
        % Setup for new level
        
        OLD_OFV                                     = LEVEL_OBJ(ixBEST);
        
        % Update options with optimal model parameters from last run
        options.POPvalues0                          = LEVEL_RESULTS(ixBEST).rawParameterInfo.fixedEffects.values;
        options.IIVvalues0                          = LEVEL_RESULTS(ixBEST).rawParameterInfo.randomEffects.values;
        options.errorParam0                         = LEVEL_RESULTS(ixBEST).rawParameterInfo.errorParameter.values;
        
        % If sequential 0/1 order absorption need to add Tlaginput 1 in between
        if absorptionModel == 3,
            POPvalues0 = NaN(1,length(options.POPvalues0)+1);
            POPvalues0([1:9 11:end]) = options.POPvalues0;
            options.POPvalues0 = POPvalues0;
            
            IIVvalues0 = NaN(1,length(options.IIVvalues0)+1);
            IIVvalues0([1:9 11:end]) = options.IIVvalues0;
            options.IIVvalues0 = IIVvalues0;
        end

        % Remove covariate from search
        covSearchBackward.covcatNamesTest(ixBEST)   = [];
        covSearchBackward.delta_forward(ixBEST)     = [];
        covSearchBackward.delta_backward(ixBEST)    = [];
        covSearchBackward.paramNamesTest(ixBEST)    = [];
        covSearchBackward.value                     = NaN(1,length(covSearchBackward.paramNamesTest));
        
        levelModel                                  = covariateModelAll{ixBEST};
        levelCount                                  = levelCount+1;
        BackwardModelName                           = ModelNameAll{ixBEST};
        
        % Update values in covSearchBackward to estimated ones in best model
        covInfo                                     = LEVEL_RESULTS(ixBEST).rawParameterInfo.covariate;
        
        for k1x = 1:length(covSearchBackward.paramNamesTest),
            paramUpdate = covSearchBackward.paramNamesTest{k1x};
            covUpdate   = covSearchBackward.covcatNamesTest{k1x};
            
            % Find parameter and covariate
            matchIX = [];
            for kkkx=1:length(covInfo.names),
                ixParam = strfind(covInfo.names{kkkx},paramUpdate);
                ixCov = strfind(covInfo.names{kkkx},covUpdate);
                if ~isempty(ixParam) && ~isempty(ixCov),
                    matchIX(end+1) = kkkx;
                end
            end
            if length(matchIX) ~= 1,
                error('Problem with getting covariate coefficient value.');
            end
            
            % Add value to covSearchBackward
            covSearchBackward.value(k1x) = covInfo.values(matchIX);
        end
    end
    
    if  isempty(covSearchBackward.paramNamesTest),
        continueBackwardSearch = 0;
        % Finalize table for this level by adding footer
        tableText(end,1:2) = {'<TF>' 'All candidate covariates removed in backward elimination.'};
    end
    
    % Export table for the level into the log file
    textDisplay = IQMconvertCellTable2ReportTable(tableText,'text');
    fprintf(fid,'%s\n',textDisplay);
end

% Get backward results
BACKWARD_OFV = OLD_OFV;
backwardModel = levelModel;

% Report backward model
fprintf(fid,'************************************\n');
fprintf(fid,'* BACKWARD MODEL RESULTS           *\n');
fprintf(fid,'************************************\n');
fprintf(fid,'\n');

fprintf(fid,'Backward model:            %s\n',BackwardModelName);
if ~isempty(backwardModel),
    fprintf(fid,'Covariates:                "%s"\n',backwardModel);
else
    fprintf(fid,'Covariates:                "NONE"\n');
end    
fprintf(fid,'Objective function value:  %g\n',BACKWARD_OFV);
fprintf(fid,'\n');

% Close report file
fclose(fid);
