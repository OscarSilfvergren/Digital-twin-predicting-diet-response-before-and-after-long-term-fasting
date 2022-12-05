LoadEverything

%% Pre-optimization
modelName = 'ModelSimpleMetabolismSteadyState';
params = ModelStartGuess(1:98);

func =@(params)ModelEstimationHealthy_Costfunction(Krssak2004_Healthy,Lerche2009_data,Magnusson1992_data,Silfvergren2021Everyday_data,time,params,modelName);

 lb = log(ParameterBounds.LowerBoundHealthy(1:98));
 ub = log(ParameterBounds.UpperBoundHealthy(1:98));

% %% Temp
% [row column] = size(ModelValidationHealthyA);
% OptimizedParamsSorted = sortrows(ModelValidationHealthyA,column);
% 
% for i = 1:row
%     optimizedParamTemp = OptimizedParamsSorted(i,1:(column-1));
%     
%     lb(1:98)   = log(optimizedParamTemp(1:98));
%     ub(1:98)   = log(optimizedParamTemp(1:98));
%     
%     lb(79:80)  = log(ParameterBounds.LowerBoundHealthy(79:80));
%     ub(79:80)  = log(ParameterBounds.UpperBoundHealthy(79:80));
%     
%     lb(85:86)  = log(ParameterBounds.LowerBoundHealthy(85:86));
%     ub(85:86)  = log(ParameterBounds.UpperBoundHealthy(85:86));
%     
%     lb(91:92)  = log(ParameterBounds.LowerBoundHealthy(91:92));
%     ub(91:92)  = log(ParameterBounds.UpperBoundHealthy(91:92));
%     
%     lb(97:98)  = log(ParameterBounds.LowerBoundHealthy(97:98));
%     ub(97:98)  = log(ParameterBounds.UpperBoundHealthy(97:98));
%     
%     
%     [ModelValidationHealthy_params, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
%     ModelValidationHealthy_params = exp(ModelValidationHealthy_params);
%     i = i
% end

%% Optimization

 [ModelValidationHealthy_params, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
 ModelValidationHealthy_params = exp(ModelValidationHealthy_params);

%% Save Best Param
save(['ModelValidationHealthy_params' datestr(now, 'yymmdd-HHMMSS')],'ModelValidationHealthy_params');

%% Kill all paralell workers and clusters

poolobj = gcp('nocreate');
myCluster = parcluster('local');
delete(myCluster.Jobs)
delete(poolobj)
