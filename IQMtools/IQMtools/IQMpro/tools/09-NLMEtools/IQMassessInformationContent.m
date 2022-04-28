function [] = IQMassessInformationContent(projectFolder,model,output,dosings,obsTimes,options)
% This function allows to predict the information content in data of (a)
% future studies, given the planned dosing and observation schedule.
%
% It is simply based on sensitivity analysis wrt to changes in the model
% parameters and correlation of the sensitivity trajectories.
%
% Assessed will be the impact of changes in single model parameters on the
% readout at given observation times. The mean of these predicted
% observations will be calculated and normalized sensitivities plotted as
% barplot. Additionally, the correlation matrix of the normalized
% sensitivities will be plotted for the parameters that have an impact on
% the readout of more than a user definable threshold.
%
% What is this function good for?
% If you have densly sampled data from Phase II studies and get a load of
% sparsly sampled Phase III data in. Then this function might support you
% in the selection of the parameters that you want to consider for
% refitting on all data together. 
%
% Covariates, random effects, residual errors are NOT taken into account!
% Population mean estimates of parameters are used!
%
% [SYNTAX]
% [] = IQMassessInformationContent(projectFolder,model,output,dosings,obsTimes)
% [] = IQMassessInformationContent(projectFolder,model,output,dosings,obsTimes,options)
%
% [INPUT]
% projectFolder:    String with the name of the NLME(NONMEM or MONOLIX) project folder for
%                   which to do the analysis. Needs to include the path to the 
% 					project folder
% model:            Structural model fitting to the NLME(NONMEM or MONOLIX) fits to use for
%                   the simulation - or path to this model
% output:           Cell-array with measured output names (model variables)
% dosings:          Cell-array with dosing schemes to simulate the model for
%                   This allows to consider more than one future study - or
%                   paths to these dosings.
% obsTimes:         Cell-array with Observation times to compare the models
%                   at. Each entry in the cell-array corresponds to the
%                   same entry in the "dosings" argument. Thus dosings and
%                   obsTimes cell-arrays should have the same length
% options:          Matlab structure with optional information
%       options.pertSize            Relative parameter perturbations are
%                                   used. (default: 10%)
%       options.sensThreshold       Sensitivity threshold to select the
%                                   parameters that have a significant
%                                   impact on the output for subsequent
%                                   correlation analysis. Default is 10%.
%                                   As example, 10% means that when
%                                   perturbing a parameter by x%, a
%                                   significant perturbation of the output
%                                   is considered to be a change of more
%                                   than x/10 percent (plus or minus).
%       options.optionsIntegrator   options for the integration.
%                                   By default: abstol=1e-6, reltol=1e-6
%       options.filename            Path and filename for export of figures
%                                   to PS (Windows) or PDF (Unix). If not
%                                   defined or empty, the figures will only
%                                   be plotted 
%
% [OUTPUT]
% Plots - which also can be exported to PDF if desired.

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

% Handle options
try optionsIntegrator           = options.optionsIntegrator;        catch, optionsIntegrator    = [];               end
try pertSize                    = options.pertSize;                 catch, pertSize             = 10;               end
try sensThreshold               = options.sensThreshold;            catch, sensThreshold        = 10;               end
try filename                    = options.filename;                 catch, filename             = '';               end
try optionsIntegrator.abstol    = optionsIntegrator.abstol;         catch, optionsIntegrator.abstol = 1e-6;         end
try optionsIntegrator.reltol    = optionsIntegrator.reltol;         catch, optionsIntegrator.reltol = 1e-6;         end

% Handle multiple dosings / obsTimes
if ischar(model),
    model = IQMmodel(model);
end
if ~iscell(dosings),
    dosings = {dosings};
end
for k=1:length(dosings),
    if ischar(dosings{k}),
        dosings{k} = IQMdosing(dosings{k});
    end
end
if ~iscell(obsTimes),
    obsTimes = {obsTimes};
end
if length(dosings) ~= length(obsTimes),
    error('Number of provided dosing scenarios and number of provided observation time vectors need to match.');
end
if ~iscell(output),
    output = {output};
end

% Generate moddos and MEX model
% Use one dosing ... need to be the same structure, only different times
% and doses are allowed
moddos      = mergemoddosIQM(model,dosings{1});
mexModel    = IQMmakeTempMEXmodel(moddos);

% Get population parameters for model
parameters = IQMsampleNLMEfitParam(projectFolder,0,0);
paramNames  = parameters.parameterNames;
paramValues = parameters.parameterValuesPopulation;

% Check that parameters in all fits are available in the model
modelparamnames = IQMparameters(moddos);
for k2=1:length(paramNames),
    ix = strmatchIQM(paramNames{k2}, modelparamnames, 'exact');
    if isempty(ix),
        error('IQMassessInformationContent: Parameters provided in the fit results ("projectFolder") need to be present in the structural model ("model").');
    end
end

% Simulate nominal PK parameters for all provided dosing schemes and 
% concatenate results 
output_nominal = [];
for k=1:length(dosings),
    % Search for Tinfinput and Tlaginput parameter names and if 0 then exchange to 1e-10
    % Otherwise numerical problems!
    ixxx                    = strmatchIQM('Tk0input',paramNames);
    for kxxx=1:length(ixxx),
        xxx = paramValues(:,ixxx(kxxx));
        xxx(xxx<=2e-10) = 1e-10;
        paramValues(:,ixxx(kxxx)) = xxx;
    end
    ixxx                    = strmatchIQM('Tinfinput',paramNames);
    for kxxx=1:length(ixxx),
        xxx = paramValues(:,ixxx(kxxx));
        xxx(xxx<=2e-10) = 1e-10;
        paramValues(:,ixxx(kxxx)) = xxx;
    end
    ix = strmatchIQM('Tlag',paramNames);
    for klag=1:length(ix),
        if paramValues(ix(klag)) < 1e-10,
            paramValues(ix(klag)) = 0;
        end
    end    
    y = IQMsimdosing(mexModel,dosings{k},obsTimes{k},[],paramNames,paramValues,optionsIntegrator);
    for k2=1:length(output),
        output_nominal_k = y.variablevalues(:,variableindexIQM(moddos,output{k2}));
        output_nominal = [output_nominal(:); output_nominal_k(:)];
    end
end

% Simulate single perturbed parameters (use all parameters in the model)
% for all provided dosing schemes and concatenate results 
output_pert    = [];
for k=1:length(paramNames),
    % Get parameter to perturb and new value
    paramName               = paramNames{k};
    paramValuePert          = paramValues(k)*(1+pertSize/100);
    % Construct new full perturbed parameter vector
    paramValues_Pert        = paramValues;
    paramValues_Pert(k)     = paramValuePert;
    
    % Simulate all dosing scenarios
    output_pert_k = [];
    for k2=1:length(dosings),
        % Search for Tinfinput and Tlaginput parameter names and if 0 then exchange to 1e-10
        % Otherwise numerical problems!
        ixxx                    = strmatchIQM('Tk0input',paramNames);
        for kxxx=1:length(ixxx),
            xxx = paramValues(:,ixxx(kxxx));
            xxx(xxx<=2e-10) = 1e-10;
            paramValues(:,ixxx(kxxx)) = xxx;
        end
        ixxx                    = strmatchIQM('Tinfinput',paramNames);
        for kxxx=1:length(ixxx),
            xxx = paramValues(:,ixxx(kxxx));
            xxx(xxx<=2e-10) = 1e-10;
            paramValues(:,ixxx(kxxx)) = xxx;
        end
        ix = strmatchIQM('Tlag',paramNames);
        for klag=1:length(ix),
            if paramValues(ix(klag)) < 1e-10,
                paramValues(ix(klag)) = 0;
            end
        end
        y                   = IQMsimdosing(mexModel,dosings{k2},obsTimes{k2},[],paramNames,paramValues_Pert,optionsIntegrator);
        for k3=1:length(output),
            output_pert_k2      = y.variablevalues(:,variableindexIQM(moddos,output{k3}));
            output_pert_k       = [output_pert_k(:); output_pert_k2(:)];
        end
    end

    % Collect output
    output_pert             = [output_pert output_pert_k];
end

% Calculate sensitivity trajectories
% Normalize by perturbation size
% Expand output_nominal
output_nominal_expanded = output_nominal(:,ones(1,length(paramNames)));
% Normalized sensitivity
Sn = 100*(output_pert - output_nominal_expanded)./output_nominal_expanded/pertSize*100;

% Calculate metric
meanSensitivity   = nanmeanIQM(abs(Sn));

% Prepare output to file if needed
IQMstartNewPrintFigure(filename);

% Display sensitivities
figure(1); clf
bar([1:length(meanSensitivity)],meanSensitivity);
paramNamesPlot = paramNames;
set(gca,'XTickLabel',paramNamesPlot)
grid on;
set(gca,'FontSize',12)
xlabel('Parameters','FontSize',14);
ylabel('Normalized sensitivities [%]','FontSize',14);
IQMprintFigure(gcf,filename)

% Assess correlation of the most important parameters
ix_important        = find(abs(meanSensitivity) >= sensThreshold);
param_corr_Names    = paramNames(ix_important);
X = [];
for k=1:length(param_corr_Names),
    ix = strmatchIQM(param_corr_Names{k},paramNames,'exact');
    X(:,k) = Sn(:,ix);
end
% Remove rows in which NaNs appear
[rowNaN,colNaN] = find(isnan(X));
X(rowNaN,:) = [];
corr_param = abs(corrcoef(X));

% Plot results
figure(2); clf;
pcolor([corr_param zeros(length(corr_param),1); zeros(1,length(corr_param)) 0])
axis square;
colorbar('EastOutside');
set(gca,'XTick',[1.5:length(param_corr_Names)+0.5]);
set(gca,'XTickLabel',param_corr_Names);
set(gca,'YTick',[1.5:length(param_corr_Names)+0.5]);
set(gca,'YTickLabel',param_corr_Names);
colormap('Bone');
set(gca,'FontSize',12)
title('Predicted correlations of parameters with significant impact','FontSize',14,'Interpreter','none');
IQMprintFigure(gcf,filename)

% Close export to file
IQMconvert2pdf(filename);
if ~isempty(filename),
    close all
end
