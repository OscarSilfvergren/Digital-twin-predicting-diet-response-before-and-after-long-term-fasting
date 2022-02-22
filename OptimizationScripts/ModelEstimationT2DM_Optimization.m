LoadEverything

%% Pre-optimization
modelName = 'ModelSimpleMetabolismSteadyState';
params = ModelStartGuess(1:98);

func =@(params)ModelEstimationT2DM_Costfunction(Krssak2004_Diabetes,Magnusson1992_data,time,params,modelName);

lb = log(ParameterBounds.LowerBoundT2D(1:98));
ub = log(ParameterBounds.UpperBoundT2D(1:98));

%% Optimization

[ModelValidationT2DM_params, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
ModelValidationT2DM_params = exp(ModelValidationT2DM_params);


%% Save Best Param
save(['ModelValidationT2DM_params' datestr(now, 'yymmdd-HHMMSS')],'ModelValidationT2DM_params');
