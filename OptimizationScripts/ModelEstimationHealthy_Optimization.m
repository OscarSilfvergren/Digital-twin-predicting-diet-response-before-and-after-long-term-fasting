LoadEverything

%% Pre-optimization
modelName = 'ModelSimpleMetabolismSteadyState';
params = ModelStartGuess(1:98);

func =@(params)ModelEstimationHealthy_Costfunction(Krssak2004_Healthy,Lerche2009_data,Magnusson1992_data,Silfvergren2021Everyday_data,time,params,modelName);

lb = log(ParameterBounds.LowerBoundHealthy(1:98));
ub = log(ParameterBounds.UpperBoundHealthy(1:98));

%% Optimization

[ModelValidationHealthy_params, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
ModelValidationHealthy_params = exp(ModelValidationHealthy_params);

%% Save Best Param
save(['ModelValidationHealthy_params' datestr(now, 'yymmdd-HHMMSS')],'ModelValidationHealthy_params');
