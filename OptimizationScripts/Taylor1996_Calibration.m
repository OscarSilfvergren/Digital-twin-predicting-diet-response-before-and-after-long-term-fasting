clear Taylor1996_ParamHealthyCalibrated

%% Pre-calc
modelName = 'ModelSimpleMetabolismSteadyState';
options2 = optimoptions(@particleswarm,'PlotFcn','pswplotbestf','MaxStallIterations',7);
params = ModelStartGuess(1:86);

body_information = [0 , 1, 180, 75.6 ];  % female, male, height, weight
meal_information = [10, 0, 138.64, 29.25];

[row column] = size(ModelValidationHealthyA);
OptimizedParamsSorted = sortrows(ModelValidationHealthyA,column);

%% Calibration

for i = 1:row
    optimizedParamTemp = OptimizedParamsSorted(i,1:(column-1));
    
    lb(1:86)   = log(optimizedParamTemp(1:86));
    ub(1:86)   = log(optimizedParamTemp(1:86));
    
    [lb,ub] = DifferentPopulationCalibration(lb,ub,ParameterBounds);
    lb(75:86)  = log(ParameterBounds.LowerBoundHealthy(75:86));
    ub(75:86)  = log(ParameterBounds.UpperBoundHealthy(75:86));
    
    func =@(params)Taylor1996_costfunction(Taylor1996_data,time,params,modelName,body_information, meal_information,2,3,3);
    [Taylor1996_ParamHealthy, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
    Taylor1996_ParamHealthy = exp(Taylor1996_ParamHealthy);
    
    %% Clean parameters; Data is of two seperate meals
    Taylor1996_ParamHealthy(87:101)    = [0,0,1,0,0,0,0,0,1,0,0,0,1,0,0];   % Steady State Meals are solid and there is no meal in t0
    Taylor1996_ParamHealthy(16) = Taylor1996_ParamHealthy(15);              % Identical meals, but consumed at different times. Params are set to not allow for population differences etc.
    Taylor1996_ParamHealthy(34) = Taylor1996_ParamHealthy(33);
    Taylor1996_ParamHealthy(56) = Taylor1996_ParamHealthy(55);
    Taylor1996_ParamHealthy(60) = Taylor1996_ParamHealthy(59);
    Taylor1996_ParamHealthy(64) = Taylor1996_ParamHealthy(63);
    Taylor1996_ParamHealthy(102:105)   = [0 , 1, 180, 75.6];                % female, male, height, weight
    Taylor1996_ParamHealthy(106:109)   = [1,0,0,0];                         % This is study A
    Taylor1996_ParamHealthy(110:117)   = [20, 0, 138.64, 29.25,1,0,0,0];    % Meal Information
    
    %% sort
    Taylor1996_ParamHealthy(column) = minCostPS;
    Taylor1996_ParamHealthyCalibrated(i,:) = Taylor1996_ParamHealthy;
    
    i = i
end

%% Save best Param
Taylor1996_ParamHealthyCalibrated = sortrows(Taylor1996_ParamHealthyCalibrated,column);
save(['Taylor1996_ParamHealthy' datestr(now, 'yymmdd-HHMMSS')],'Taylor1996_ParamHealthyCalibrated');
