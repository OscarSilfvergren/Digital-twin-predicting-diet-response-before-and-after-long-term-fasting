
%% Pre-calc
modelName = 'ModelSimpleMetabolismFast';
model = IQMmodel([modelName '.txt']);
IQMmakeMEXmodel(model);
model = str2func(modelName);
options2 = optimoptions(@particleswarm,'PlotFcn','pswplotbestf','MaxStallIterations',7);

time = linspace(0,11500,11501);
params = ModelStartGuess(1:80);

[row column] = size(ModelValidationHealthyA);
OptimizedParamsSorted = sortrows(ModelValidationHealthyA,column);

%% Fasting intervention

clear Silfvergren2021_ParamHealthyCalibratedp1
clear Silfvergren2021_ParamHealthyCalibratedp2

for i = 1:row
    
    optimizedParamTemp = OptimizedParamsSorted(i,1:(column-1));
    
    lb(1:80)   = log(optimizedParamTemp(1:80));
    ub(1:80)   = log(optimizedParamTemp(1:80));

    [lb,ub] = DifferentPopulationCalibration(lb,ub,ParameterBounds);
    lb(75:80)  = log(ParameterBounds.LowerBoundHealthy(75:80));
    ub(75:80)  = log(ParameterBounds.UpperBoundHealthy(75:80));
    
    % P1
    body_information = [1 , 0, 178, 84 ];  % female, male, height, weight
    meal_information = [1, 0, 0, 0];
    
    func =@(params)Silfvergren2021_costfunction(Silfvergren2021_data.p1_glucoseCalibrated,Silfvergren2021_data.time,time,params,modelName,body_information,meal_information,3);
    [Silfvergren2021_ParamHealthy, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
    Silfvergren2021_ParamHealthy = exp(Silfvergren2021_ParamHealthy);
    Silfvergren2021_ParamHealthy = AssignParameter(Silfvergren2021_ParamHealthy, body_information, meal_information);
    Silfvergren2021_ParamHealthy(column) = minCostPS;
    Silfvergren2021_ParamHealthyCalibratedp1(i,:) = Silfvergren2021_ParamHealthy;
    
    % P2
    body_information = [0 , 1, 180, 87 ];  % female, male, height, weight
    meal_information = [1, 0, 0, 0];
    
    func =@(params)Silfvergren2021_costfunction(Silfvergren2021_data.p2_glucose,Silfvergren2021_data.time,time,params,modelName,body_information,meal_information,4);
    [Silfvergren2021_ParamHealthy, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
    Silfvergren2021_ParamHealthy = exp(Silfvergren2021_ParamHealthy);
    Silfvergren2021_ParamHealthy = AssignParameter(Silfvergren2021_ParamHealthy, body_information, meal_information);
    Silfvergren2021_ParamHealthy(column) = minCostPS;
    Silfvergren2021_ParamHealthyCalibratedp2(i,:) = Silfvergren2021_ParamHealthy;
    
    i = i
end

% Save best Param
Silfvergren2021_ParamHealthyCalibratedp1 = sortrows(Silfvergren2021_ParamHealthyCalibratedp1,column);
save(['Silfvergren2021_ParamHealthyCalibratedp1' datestr(now, 'yymmdd-HHMMSS')],'Silfvergren2021_ParamHealthyCalibratedp1');

Silfvergren2021_ParamHealthyCalibratedp2 = sortrows(Silfvergren2021_ParamHealthyCalibratedp2,column);
save(['Silfvergren2021_ParamHealthyCalibratedp2' datestr(now, 'yymmdd-HHMMSS')],'Silfvergren2021_ParamHealthyCalibratedp2');

%% OPTT; P1

clear Silfvergren2021_ParamP1fedCalibrated
clear Silfvergren2021_ParamP1unfedCalibrated

for i = 1:row
    optimizedParamTemp = OptimizedParamsSorted(i,1:(column-1));
    
    lb(1:80)   = log(optimizedParamTemp(1:80));
    ub(1:80)   = log(optimizedParamTemp(1:80));

    [lb,ub] = DifferentPopulationCalibration(lb,ub,ParameterBounds);
    lb(75:80)  = log(ParameterBounds.LowerBoundHealthy(75:80));
    ub(75:80)  = log(ParameterBounds.UpperBoundHealthy(75:80));
    
    body_information = [1 , 0, 178, 84 ];  % female, male, height, weight
    meal_information = [3, 0, 2.6, 25.55];
    
    % fed
    func =@(params)Silfvergren2021_costfunction(Silfvergren2021_data.value_fed_p1(1:24),Silfvergren2021_data.time_fed_p1(1:24),time,params,modelName,body_information,meal_information,2);
    [Silfvergren2021_ParamP1fed, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
    Silfvergren2021_ParamP1fed = exp(Silfvergren2021_ParamP1fed);
    Silfvergren2021_ParamP1fed = AssignParameter(Silfvergren2021_ParamP1fed, body_information, meal_information);
    Silfvergren2021_ParamP1fed(column) = minCostPS;
    Silfvergren2021_ParamP1fedCalibrated(i,:) = Silfvergren2021_ParamP1fed;
    
    % unfed
    func =@(params)Silfvergren2021_costfunction(Silfvergren2021_data.Value_fasted_p1(1:60),Silfvergren2021_data.time_fastedStart_p1(1:59),time,params,modelName,body_information,meal_information,2);
    [Silfvergren2021_ParamP1unfed, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
    Silfvergren2021_ParamP1unfed = exp(Silfvergren2021_ParamP1unfed);
    Silfvergren2021_ParamP1unfed = AssignParameter(Silfvergren2021_ParamP1unfed, body_information, meal_information);
    Silfvergren2021_ParamP1unfed(column) = minCostPS;
    Silfvergren2021_ParamP1unfedCalibrated(i,:) = Silfvergren2021_ParamP1unfed;
    
    i = i
end

% Save best Param
Silfvergren2021_ParamP1fedCalibrated = sortrows(Silfvergren2021_ParamP1fedCalibrated,column);
save(['Silfvergren2021_ParamP1fedCalibrated' datestr(now, 'yymmdd-HHMMSS')],'Silfvergren2021_ParamP1fedCalibrated');

Silfvergren2021_ParamP1unfedCalibrated = sortrows(Silfvergren2021_ParamP1unfedCalibrated,column);
save(['Silfvergren2021_ParamP1unfedCalibrated' datestr(now, 'yymmdd-HHMMSS')],'Silfvergren2021_ParamP1unfedCalibrated');

%% OPTT; P2

clear Silfvergren2021_ParamP2fedCalibrated
clear Silfvergren2021_ParamP2unfedCalibrated

for i = 1:row
    optimizedParamTemp = OptimizedParamsSorted(i,1:(column-1));
    lb(1:80)   = log(optimizedParamTemp(1:80));
    ub(1:80)   = log(optimizedParamTemp(1:80));

    [lb,ub] = DifferentPopulationCalibration(lb,ub,ParameterBounds);
    lb(75:80)  = log(ParameterBounds.LowerBoundHealthy(75:80));
    ub(75:80)  = log(ParameterBounds.UpperBoundHealthy(75:80));
    
    body_information = [1 , 0, 180, 87 ];  % female, male, height, weight
    meal_information = [3, 0, 2.6, 25.55];
    
    % fed
    func =@(params)Silfvergren2021_costfunction(Silfvergren2021_data.value_fed_p2(1:21),Silfvergren2021_data.time_fed_p2(1:21),time,params,modelName,body_information,meal_information,2);
    [Silfvergren2021_ParamP2fed, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
    Silfvergren2021_ParamP2fed = exp(Silfvergren2021_ParamP2fed);
    Silfvergren2021_ParamP2fed = AssignParameter(Silfvergren2021_ParamP2fed, body_information, meal_information);
    Silfvergren2021_ParamP2fed(column) = minCostPS;
    Silfvergren2021_ParamP2fedCalibrated(i,:) = Silfvergren2021_ParamP2fed;
    
    % unfed
    func =@(params)Silfvergren2021_costfunction(Silfvergren2021_data.Value_fasted_p2(1:46),Silfvergren2021_data.time_fastedStart_p2(1:44),time,params,modelName,body_information,meal_information,3);
    [Silfvergren2021_ParamP2unfed, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
    Silfvergren2021_ParamP2unfed = exp(Silfvergren2021_ParamP2unfed);
    Silfvergren2021_ParamP2unfed = AssignParameter(Silfvergren2021_ParamP2unfed, body_information, meal_information);
    Silfvergren2021_ParamP2unfed(column) = minCostPS;
    Silfvergren2021_ParamP2unfedCalibrated(i,:) = Silfvergren2021_ParamP2unfed;
    
    i = i
end

% Save best Param
Silfvergren2021_ParamP2fedCalibrated = sortrows(Silfvergren2021_ParamP2fedCalibrated,column);
save(['Silfvergren2021_ParamP2fedCalibrated' datestr(now, 'yymmdd-HHMMSS')],'Silfvergren2021_ParamP2fedCalibrated');

Silfvergren2021_ParamP2unfedCalibrated = sortrows(Silfvergren2021_ParamP2unfedCalibrated,column);
save(['Silfvergren2021_ParamP2unfedCalibrated' datestr(now, 'yymmdd-HHMMSS')],'Silfvergren2021_ParamP2unfedCalibrated');
