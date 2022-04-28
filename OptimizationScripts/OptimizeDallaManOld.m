LoadEverything

modelName = 'DallamanOldData';
model = IQMmodel([modelName '.txt']);
IQMmakeMEXmodel(model);
model = str2func(modelName);
[pNames, param] = IQMparameters(model);

%% Pre-optimization

params = param(1:40);

func = @(params)OldModelDallaMan2007_costfunctionHealthy(DallaMan2007_data,params,modelName);

lb = log(param(1:40)*0.8);
ub = log(param(1:40)*1.2);

%% Optimization

[OptimizeDallaManOld_params, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
OptimizeDallaManOld_params = exp(OptimizeDallaManOld_params);

%% Save Best Param
save(['OptimizeDallaManOld_params' datestr(now, 'yymmdd-HHMMSS')],'OptimizeDallaManOld_params');
