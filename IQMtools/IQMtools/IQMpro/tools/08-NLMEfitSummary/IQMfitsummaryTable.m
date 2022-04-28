function [tableCell] = IQMfitsummaryTable(projectPaths,filename)
% This function creates a typical report-type NLME model parameter table.
% "projectPaths" defines the NLME projects that are to be reported in this
% table.
%
% [SYNTAX]
% [tableCell] = IQMfitsummaryTable(projectPaths)
% [tableCell] = IQMfitsummaryTable(projectPaths,filename)
%
% [INPUT]
% projectPaths:     Cell-array with paths to NLME projects to create the
%                   fit results table for. The table reports each model
%                   result in a column. Columns ordered according to the
%                   order in projectPaths.
% filename:         Path and filename where to generate the output file.
%                   default: '' (no file generated).
%
% [OUTPUT]
% Table shown in command window and if desired also saved to file for
% reporting purpose.

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

% Handle variable input arguments
if nargin<2,
    filename = '';
end

% Handle cell-array thing
if ischar(projectPaths),
    projectPaths = {projectPaths};
end

% Check projects if NLME
for k=1:length(projectPaths),
    if ~isNLMEprojectIQM(projectPaths{k}),
        error('Path "%s" does not contain an NLME project.',projectPaths{k});
    end
end

% Load project results - do not reorder
order   = '';
RESULTS = parseSelectedProjectFolderResultsIQM(projectPaths,order);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process the RESULT information somewhat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine all available fixed effect parameters in the models
ALLfixEffectNames = {};
for k=1:length(RESULTS),
    if ~isempty(RESULTS(k).rawParameterInfo),
        fek = RESULTS(k).rawParameterInfo.fixedEffects.names;
        ALLfixEffectNames = [ALLfixEffectNames fek];
    end
end
ALLfixEffectNames = unique(ALLfixEffectNames);

% Determine all available random effect parameters in the models
ALLrandomEffectNames = {};
for k=1:length(RESULTS),
    if ~isempty(RESULTS(k).rawParameterInfo),
        fek = RESULTS(k).rawParameterInfo.randomEffects.names;
        ALLrandomEffectNames = [ALLrandomEffectNames fek];
    end
end
ALLrandomEffectNames = unique(ALLrandomEffectNames);

% Determine all available correlation parameters in the models
ALLcorrelationNames = {};
for k=1:length(RESULTS),
    if ~isempty(RESULTS(k).rawParameterInfo),
        fek = RESULTS(k).rawParameterInfo.correlation.names;
        ALLcorrelationNames = [ALLcorrelationNames fek];
    end
end
ALLcorrelationNames = unique(ALLcorrelationNames);

% Determine all available covariate parameters in the models
ALLcovariateNames = {};
for k=1:length(RESULTS),
    if ~isempty(RESULTS(k).rawParameterInfo),
        fek = RESULTS(k).rawParameterInfo.covariate.names;
        ALLcovariateNames = [ALLcovariateNames fek];
    end
end
ALLcovariateNames = unique(ALLcovariateNames);

% Determine all available covariate parameters in the models
ALLerrorNames = {};
for k=1:length(RESULTS),
    if ~isempty(RESULTS(k).rawParameterInfo),
        fek = RESULTS(k).rawParameterInfo.errorParameter.names;
        ALLerrorNames = [ALLerrorNames fek];
    end
end
ALLerrorNames = unique(ALLerrorNames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Retain last name of model only
modelsNames_Short = {RESULTS.model};
for k=1:length(modelsNames_Short),
    [~,f,e] = fileparts(modelsNames_Short{k});
    modelsNames_Short{k} = [f e];
end
    


tableCell                               = {'<TT>' 'Population parameter estimates for considered models'};
tableCell(end+1,1:2+length(RESULTS))    = {'<TH>' ''          modelsNames_Short{:}};
tableCell(end+1,1:2)                    = {'<TR>' 'Parameter'};
for k=1:length(RESULTS),
    tableCell(end,2+k) = {'Value      (RSE%)        (trans)'};
end

% FIXED EFFECTs - non-estimated get postfix (FIX).
FLAG_FOOTER = 0;
for k0=1:length(ALLfixEffectNames),
    tableCell(end+1,1:2) = {'<TR>' ALLfixEffectNames{k0}};
    for k=1:length(RESULTS),
        % Check if parameter in model
        ix = strmatchIQM(ALLfixEffectNames{k0},RESULTS(k).rawParameterInfo.fixedEffects.names,'exact');
        if ~isempty(ix),
            modelVALUES             = RESULTS(k).rawParameterInfo.fixedEffects.values(ix);
            modelRSES               = RESULTS(k).rawParameterInfo.fixedEffects.rse(ix);
            modelESTIMATED          = RESULTS(k).rawParameterInfo.fixedEffects.estimated(ix);
            if modelESTIMATED,
                value               = sprintf('%1.3g',modelVALUES);
                if RESULTS(k).NONMEM,
                    rse             = sprintf('(%1.3g%%*)',modelRSES);
                    FLAG_FOOTER     = 1;
                else
                    rse                 = sprintf('(%1.3g%%)',modelRSES);
                end
                valuerse            = sprintf('%s %s',postFillCharIQM(value,10,' '),postFillCharIQM(rse,14,' '));
            else
                if ~isnan(modelVALUES),
                    value           = sprintf('%1.3g',modelVALUES);
                    rse             = sprintf('(FIX)');
                    valuerse        = sprintf('%s %s',postFillCharIQM(value,10,' '),postFillCharIQM(rse,14,' '));
                else
                    valuerse        = '-';
                end
            end
            % Get MU referencing transformation of fixed effect
            modelTRANS              = RESULTS(k).rawParameterInfo.fixedEffects.distribution_info{ix};
            if strcmp(modelTRANS,'(psi)'),
                modelTRANS = '(norm)';
            elseif strcmp(modelTRANS,'log(psi./(1-psi))'),
                modelTRANS = '(logit)';
            elseif strcmp(modelTRANS,'log(psi)'),
                modelTRANS = '(log)';
            else
                error('Unknown transformation of fixed effect.');
            end
            % Add to text
            valuerse = sprintf('%s%s',valuerse,modelTRANS);
           
        else
            valuerse = '-';
        end
        % Add to table
        tableCell(end,2+k) = {valuerse};
    end
end

% Add footer if needed
if FLAG_FOOTER,
    tableCell(end+1,1:2) = {'<TF>' '* Standard errors for NONMEM models approximated by sampling due to MU-Referencing'};    
end

% Separator
tableCell(end+1,1) = {'<HR>'};

% RANDOM EFFECTs - non-estimated get postfix (FIX).
for k0=1:length(ALLrandomEffectNames),
    tableCell(end+1,1:2) = {'<TR>' ALLrandomEffectNames{k0}};
    for k=1:length(RESULTS),
        % Check if parameter in model
        ix = strmatchIQM(ALLrandomEffectNames{k0},RESULTS(k).rawParameterInfo.randomEffects.names,'exact');
        if ~isempty(ix),
            modelVALUES             = RESULTS(k).rawParameterInfo.randomEffects.values(ix);
            modelRSES               = RESULTS(k).rawParameterInfo.randomEffects.rse(ix);
            modelESTIMATED          = RESULTS(k).rawParameterInfo.randomEffects.estimated(ix);
            if modelESTIMATED==1,
                value               = sprintf('%1.3g',modelVALUES);
                rse                 = sprintf('(%1.3g%%)',modelRSES);
                valuerse            = sprintf('%s %s',postFillCharIQM(value,10,' '),rse);
            else
                if ~isnan(modelVALUES),
                    value           = sprintf('%1.3g',modelVALUES);
                    rse             = sprintf('(FIX)',modelRSES);
                    valuerse        = sprintf('%s %s',postFillCharIQM(value,10,' '),rse);
                else
                    valuerse        = '-';
                end
            end
        else
            valuerse = '-';
        end
        % Add to table
        tableCell(end,2+k) = {valuerse};
    end
end
      
% Separator
tableCell(end+1,1) = {'<HR>'};

% CORRELATIONS - non-estimated get postfix (FIX).
for k0=1:length(ALLcorrelationNames),
    tableCell(end+1,1:2) = {'<TR>' ALLcorrelationNames{k0}};
    for k=1:length(RESULTS),
        % Check if parameter in model
        ix = strmatchIQM(ALLcorrelationNames{k0},RESULTS(k).rawParameterInfo.correlation.names,'exact');
        if ~isempty(ix),
            modelVALUES             = RESULTS(k).rawParameterInfo.correlation.values(ix);
            modelRSES               = RESULTS(k).rawParameterInfo.correlation.rse(ix);
            modelESTIMATED          = RESULTS(k).rawParameterInfo.correlation.estimated(ix);
            if modelESTIMATED,
                value               = sprintf('%1.3g',modelVALUES);
                rse                 = sprintf('(%1.3g%%)',modelRSES);
                valuerse            = sprintf('%s %s',postFillCharIQM(value,10,' '),rse);
            else
                if ~isnan(modelVALUES),
                    value           = sprintf('%1.3g',modelVALUES);
                    rse             = sprintf('(FIX)',modelRSES);
                    valuerse        = sprintf('%s %s',postFillCharIQM(value,10,' '),rse);
                else
                    valuerse        = '-';
                end
            end
        else
            valuerse = '-';
        end
        % Add to table
        tableCell(end,2+k) = {valuerse};
    end
end
      
% Separator
tableCell(end+1,1) = {'<HR>'};

% COVARIATES - non-estimated get postfix (FIX).
for k0=1:length(ALLcovariateNames),
    tableCell(end+1,1:2) = {'<TR>' ALLcovariateNames{k0}};
    for k=1:length(RESULTS),
        % Check if parameter in model
        ix = strmatchIQM(ALLcovariateNames{k0},RESULTS(k).rawParameterInfo.covariate.names,'exact');
        if ~isempty(ix),
            modelVALUES             = RESULTS(k).rawParameterInfo.covariate.values(ix);
            modelRSES               = RESULTS(k).rawParameterInfo.covariate.rse(ix);
            modelESTIMATED          = RESULTS(k).rawParameterInfo.covariate.estimated(ix);
            if modelESTIMATED,
                value               = sprintf('%1.3g',modelVALUES);
                rse                 = sprintf('(%1.3g%%)',modelRSES);
                valuerse            = sprintf('%s %s',postFillCharIQM(value,10,' '),postFillCharIQM(rse,14,' '));
            else
                if ~isnan(modelVALUES),
                    value           = sprintf('%1.3g',modelVALUES);
                    rse             = sprintf('(FIX)',modelRSES);
                    valuerse        = sprintf('%s %s',postFillCharIQM(value,10,' '),postFillCharIQM(rse,14,' '));
                else
                    valuerse        = '-';
                end
            end
            transinfo = getCovCatTransInfo(ALLcovariateNames{k0},RESULTS(k));
            valuerse = sprintf('%s%s',valuerse,transinfo);
        else
            valuerse = '-';
        end
        % Add to table
        tableCell(end,2+k) = {valuerse};
    end
end
      
% Separator
tableCell(end+1,1) = {'<HR>'};

% ERROR PARAMETERS - non-estimated get postfix (FIX).
for k0=1:length(ALLerrorNames),
    tableCell(end+1,1:2) = {'<TR>' ALLerrorNames{k0}};
    for k=1:length(RESULTS),
        % Check if parameter in model
        ix = strmatchIQM(ALLerrorNames{k0},RESULTS(k).rawParameterInfo.errorParameter.names,'exact');
        if ~isempty(ix),
            modelVALUES             = RESULTS(k).rawParameterInfo.errorParameter.values(ix);
            modelRSES               = RESULTS(k).rawParameterInfo.errorParameter.rse(ix);
            modelESTIMATED          = RESULTS(k).rawParameterInfo.errorParameter.estimated(ix);
            if modelESTIMATED,
                value               = sprintf('%1.3g',modelVALUES);
                rse                 = sprintf('(%1.3g%%)',modelRSES);
                valuerse            = sprintf('%s %s',postFillCharIQM(value,10,' '),rse);
            else
                if ~isnan(modelVALUES),
                    value           = sprintf('%1.3g',modelVALUES);
                    rse             = sprintf('(FIX)',modelRSES);
                    valuerse        = sprintf('%s %s',postFillCharIQM(value,10,' '),rse);
                else
                    valuerse        = '-';
                end
            end
        else
            valuerse = '-';
        end
        % Add to table
        tableCell(end,2+k) = {valuerse};
    end
end
      
% Separator
tableCell(end+1,1) = {'<HR>'};

% OBJ, AIC, BIC
tableCell(end+1,:) = {'<TR>' 'OBJ' RESULTS.OBJ};
tableCell(end+1,:) = {'<TR>' 'AIC' RESULTS.AIC};
tableCell(end+1,:) = {'<TR>' 'BIC' RESULTS.BIC};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display results and save to file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert to text and display text 
textDisplay = IQMconvertCellTable2ReportTable(tableCell,'text');
disp(textDisplay);

% Convert to report text and export to file if filename defined
IQMconvertCellTable2ReportTable(tableCell,'report',filename);     



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to generate transformation and reference value information for
% covariates ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [transinfo] = getCovCatTransInfo(ALLcovariateNamesK0,RESULTS)
COVNAMES            = RESULTS.projectHeader.COVNAMES;
CATNAMES            = RESULTS.projectHeader.CATNAMES;
BETACOVNAMES        = RESULTS.projectHeader.BETACOVNAMES;
BETACOVTRANS        = RESULTS.projectHeader.BETACOVTRANS;
BETACATNAMES        = RESULTS.projectHeader.BETACATNAMES;
BETACATREFERENCE    = RESULTS.projectHeader.BETACATREFERENCE;

% Process covariate name to get parameter and covariate name
name    = strrep(ALLcovariateNamesK0,'beta_','');
name    = strrep(name,')','');
name    = strrep(name,'(',',');
terms   = explodePCIQM(name);
pn      = terms{1};
cn      = terms{2};
% Check if categorical or continuous
ctypecont = NaN;
ix = strmatchIQM(cn,COVNAMES,'exact');
if ~isempty(ix),
    ctypecont = 1;
    cc = '';
else
    terms   = explodePCIQM(cn,'_');
    cn = terms{1};
    cc = terms{2};
    ix = strmatchIQM(cn,CATNAMES,'exact');
    if ~isempty(ix),
        ctypecont = 0;
    end
end
if isnan(ctypecont),
    error('Difficulty to decide covariate type.');
end
% Handle continuous covariate
if ctypecont,
    % Only need to get the covariate transformation
    ix = strmatchIQM(ALLcovariateNamesK0,BETACOVNAMES,'exact');
    trans = BETACOVTRANS{ix};
    trans = strrep(trans,'cov',cn);
%     transinfo = strrep(sprintf('(%s)',trans),'.','');
    transinfo = sprintf('(%s)',trans);
else
    % Need to get reference category
    test = sprintf('beta_%s(%s)',pn,cn);
    ix = strmatchIQM(test,BETACATNAMES,'exact');
    trans = BETACATREFERENCE{ix};
    transinfo = sprintf('(Ref: %s)',trans);
end
return
