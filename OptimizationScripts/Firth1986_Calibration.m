clear Firth1986_ParamHealthyCalibrated

%% Pre-calculation
modelName = 'ModelSimpleMetabolismSteadyState';

options2 = optimoptions(@particleswarm,'PlotFcn','pswplotbestf','MaxStallIterations',7);
params = ModelStartGuess(1:80);

body_information = [1 , 0, 175, 90 ];  % female, male, height, weight
meal_information = [10, 0, 50, 0];

[row column] = size(ModelValidationHealthyA);
OptimizedParamsSorted = sortrows(ModelValidationHealthyA,column);

for i = 1:row
    optimizedParamTemp = OptimizedParamsSorted(i,1:(column-1));
    
    lb(1:80)   = log(optimizedParamTemp(1:80));
    ub(1:80)   = log(optimizedParamTemp(1:80));

    [lb,ub] = DifferentPopulationCalibration(lb,ub,ParameterBounds);
    lb(75:80)  = log(ParameterBounds.LowerBoundHealthy(75:80));
    ub(75:80)  = log(ParameterBounds.UpperBoundHealthy(75:80));
    
    %% Sim
    func =@(params)Firth1986_costfunctionHealthy(Firth1986_data,time,params,modelName,meal_information,body_information,3,5,5);
    
    [Firth1986_ParamHealthy, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
    Firth1986_ParamHealthy = exp(Firth1986_ParamHealthy);
    Firth1986_ParamHealthy = AssignParameter(Firth1986_ParamHealthy, body_information, meal_information);
    Firth1986_ParamHealthy(column) = minCostPS;
    Firth1986_ParamHealthyCalibrated(i,:) = Firth1986_ParamHealthy;
    i = i
end

Firth1986_ParamHealthyCalibrated = sortrows(Firth1986_ParamHealthyCalibrated,column);
save(['Firth1986_ParamHealthyCalibrated' datestr(now, 'yymmdd-HHMMSS')],'Firth1986_ParamHealthyCalibrated');
