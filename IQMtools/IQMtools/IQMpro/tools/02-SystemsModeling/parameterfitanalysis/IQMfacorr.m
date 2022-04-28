function [] = IQMfacorr(estdata, varargin)
% IQMfacorr: This function determines the correlation matrix for the
% parameter sets determined with the IQMparameterfitanalysis function.
% The closer the magnitude of the values is to one, the more correlated the
% parameters. 
%
% Results are generated only for the global parameters.
%
% USAGE:
% ======
% IQMfacorr(estdata)
%
% estdata:  The estimation data returned by the function IQMparameterfitanalysis

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

% Get the parameter information
parameters = estdata.parameters;
G = estdata.Popt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DETERMINE THE PARAMETER CORRELATION MATRIX
% Take out parameters with zero variance!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n,m] = size(G);
C = cov(G);
zerovarianceindices = find(diag(C)==0);
G(:,zerovarianceindices) = [];  % take out these parameters
allparameters = parameters;
parameters = parameters(setdiff([1:length(parameters)],zerovarianceindices));
[correlationMatrix,P,LB,UB] = corrcoef(G);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY NOTE IF PARAMTERS HAVE BEEN TAKEN OUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(zerovarianceindices),
    text = '';
    for k=1:length(zerovarianceindices),
        text = sprintf('%sParameter ''%s'' shows 0 variance. Taken out of consideration.\n',text,allparameters{zerovarianceindices(k)});
    end
    disp(text);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Handle the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the correlation matrix (absolute values)
% Prepare plot matrix
plotMatrix = [correlationMatrix zeros(size(correlationMatrix,1),1); 0*ones(1,size(correlationMatrix,2)+1)];
plotMatrix = abs(plotMatrix);
% Plot the result
figH = figure; clf;
axesH = gca(figH);
pcolor(plotMatrix);
axis square;
colorbar('EastOutside','YTick',[-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1]);
set(axesH,'XTick',[1.5:size(correlationMatrix,1)+0.5]);
set(axesH,'XTickLabel',parameters);
try
    set(gca,'XTickLabelRotation',45);
catch
end
set(axesH,'YTick',[1.5:size(correlationMatrix,1)+0.5]);
set(axesH,'YTickLabel',parameters);
colormap('Bone');
title('Parameter Correlation Matrix (absolute values)');
