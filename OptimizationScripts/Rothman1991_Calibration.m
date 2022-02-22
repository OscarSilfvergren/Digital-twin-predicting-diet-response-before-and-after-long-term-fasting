
%% Pre-calc
modelName = 'ModelSimpleMetabolismSteadyState';
options2 = optimoptions(@particleswarm,'PlotFcn','pswplotbestf','MaxStallIterations',7);
params = ModelStartGuess(1:80);

body_information = [0 , 1, 180, 75 ];  % female, male, height, weight
meal_information = [1, 0, 0, 0];

[row column] = size(ModelValidationHealthyA);
OptimizedParamsSorted = sortrows(ModelValidationHealthyA,column);

for i = 1:row
    %% Bounds
    optimizedParamTemp = OptimizedParamsSorted(i,1:(column-1));
    
    lb(1:74)   = log(optimizedParamTemp(1:74));
    ub(1:74)   = log(optimizedParamTemp(1:74));
    
    [lb,ub] = DifferentPopulationCalibration(lb,ub,ParameterBounds); 
    lb(75:80)  = log(ParameterBounds.LowerBoundHealthy(75:80));
    ub(75:80)  = log(ParameterBounds.UpperBoundHealthy(75:80));
    
    %% p1
    func =@(params)Rothman1991_costfunctionIndividual(Rothman1991_data.glycogenp1(1:1),Rothman1991_data.timep1(1:1),time,params,modelName,body_information, meal_information);
    [Rothman1991_paramP1, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
    Rothman1991_paramP1 = exp(Rothman1991_paramP1);
    Rothman1991_paramP1 = AssignParameter(Rothman1991_paramP1, body_information, meal_information);
    Rothman1991_paramP1(column) = minCostPS;
    Rothman1991_paramP1Calibrated(i,:) = Rothman1991_paramP1;
    
    %% P2
    func =@(params)Rothman1991_costfunctionIndividual(Rothman1991_data.glycogenp2(1:1),Rothman1991_data.timep2(1:1),time,params,modelName,body_information, meal_information);
    [Rothman1991_paramP2, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
    Rothman1991_paramP2 = exp(Rothman1991_paramP2);
    Rothman1991_paramP2 = AssignParameter(Rothman1991_paramP2, body_information, meal_information);
    Rothman1991_paramP2(column) = minCostPS;
    Rothman1991_paramP2Calibrated(i,:) = Rothman1991_paramP2;
    
    %% P3
    func =@(params)Rothman1991_costfunctionIndividual(Rothman1991_data.glycogenp3(1:1),Rothman1991_data.timep3(1:1),time,params,modelName,body_information, meal_information);
    [Rothman1991_paramP3, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
    Rothman1991_paramP3 = exp(Rothman1991_paramP3);
    Rothman1991_paramP3 = AssignParameter(Rothman1991_paramP3, body_information, meal_information);
    Rothman1991_paramP3(column) = minCostPS;
    Rothman1991_paramP3Calibrated(i,:) = Rothman1991_paramP3;
    
    %% P4
    func =@(params)Rothman1991_costfunctionIndividual(Rothman1991_data.glycogenp4(1:1),Rothman1991_data.timep4(1:1),time,params,modelName,body_information, meal_information);
    [Rothman1991_paramP4, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
    Rothman1991_paramP4 = exp(Rothman1991_paramP4);
    Rothman1991_paramP4 = AssignParameter(Rothman1991_paramP4, body_information, meal_information);
    Rothman1991_paramP4(column) = minCostPS;
    Rothman1991_paramP4Calibrated(i,:) = Rothman1991_paramP4;
    
    %% P5
    func =@(params)Rothman1991_costfunctionIndividual(Rothman1991_data.glycogenp5(1:1),Rothman1991_data.timep5(1:1),time,params,modelName,body_information, meal_information);
    [Rothman1991_paramP5, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
    Rothman1991_paramP5 = exp(Rothman1991_paramP5);
    Rothman1991_paramP5 = AssignParameter(Rothman1991_paramP5, body_information, meal_information);
    Rothman1991_paramP5(column) = minCostPS;
    Rothman1991_paramP5Calibrated(i,:) = Rothman1991_paramP5;
    
    %% P6
    func =@(params)Rothman1991_costfunctionIndividual(Rothman1991_data.glycogenp6(1:1),Rothman1991_data.timep6(1:1),time,params,modelName,body_information, meal_information);
    [Rothman1991_paramP6, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
    Rothman1991_paramP6 = exp(Rothman1991_paramP6);
    Rothman1991_paramP6 = AssignParameter(Rothman1991_paramP6, body_information, meal_information);
    Rothman1991_paramP6(column) = minCostPS;
    Rothman1991_paramP6Calibrated(i,:) = Rothman1991_paramP6;
    
    %% P7
    func =@(params)Rothman1991_costfunctionIndividual(Rothman1991_data.glycogenp7(1:1),Rothman1991_data.timep7(1:1),time,params,modelName,body_information, meal_information);
    [Rothman1991_paramP7, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
    Rothman1991_paramP7 = exp(Rothman1991_paramP7);
    Rothman1991_paramP7 = AssignParameter(Rothman1991_paramP7, body_information, meal_information);
    Rothman1991_paramP7(column) = minCostPS;
    Rothman1991_paramP7Calibrated(i,:) = Rothman1991_paramP7;
    
    %% Mean
    func =@(params)Rothman1991_costfunctionMean(Rothman1991_data,time,params,modelName,body_information, meal_information,1,1);
    [Rothman1991_paramMean, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
    Rothman1991_paramMean = exp(Rothman1991_paramMean);
    Rothman1991_paramMean = AssignParameter(Rothman1991_paramMean, body_information, meal_information);
    Rothman1991_paramMean(column) = minCostPS;
    Rothman1991_paramPMeanCalibrated(i,:) = Rothman1991_paramMean;
    
    i = i
end

%% Save
Rothman1991_paramP1Calibrated = sortrows(Rothman1991_paramP1Calibrated,column);
save(['Rothman1991P1' datestr(now, 'yymmdd-HHMMSS')],'Rothman1991_paramP1Calibrated');

Rothman1991_paramP2Calibrated = sortrows(Rothman1991_paramP2Calibrated,column);
save(['Rothman1991P2' datestr(now, 'yymmdd-HHMMSS')],'Rothman1991_paramP2Calibrated');

Rothman1991_paramP3Calibrated = sortrows(Rothman1991_paramP3Calibrated,column);
save(['Rothman1991P3' datestr(now, 'yymmdd-HHMMSS')],'Rothman1991_paramP3Calibrated');

Rothman1991_paramP4Calibrated = sortrows(Rothman1991_paramP4Calibrated,column);
save(['Rothman1991P4' datestr(now, 'yymmdd-HHMMSS')],'Rothman1991_paramP4Calibrated');

Rothman1991_paramP5Calibrated = sortrows(Rothman1991_paramP5Calibrated,column);
save(['Rothman1991P5' datestr(now, 'yymmdd-HHMMSS')],'Rothman1991_paramP5Calibrated');

Rothman1991_paramP6Calibrated = sortrows(Rothman1991_paramP6Calibrated,column);
save(['Rothman1991P6' datestr(now, 'yymmdd-HHMMSS')],'Rothman1991_paramP6Calibrated');

Rothman1991_paramP7Calibrated = sortrows(Rothman1991_paramP7Calibrated,column);
save(['Rothman1991P7' datestr(now, 'yymmdd-HHMMSS')],'Rothman1991_paramP7Calibrated');

Rothman1991_paramPMeanCalibrated = sortrows(Rothman1991_paramPMeanCalibrated,column);
save(['Rothman1991Mean' datestr(now, 'yymmdd-HHMMSS')],'Rothman1991_paramPMeanCalibrated');


