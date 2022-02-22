function Cost_Model = Taylor1996_costfunction(Taylor1996_data,time,params,modelName,body_information, meal_information,datapointsGly,datapointsG,datapointsI)

model = str2func(modelName);
params=exp(params);

params(16) = params(15); % Identical meals, but consumed at different times. Params are set to not allow differences in insulin response etc.
params(34) = params(33);
params(56) = params(55);
params(60) = params(59);
params(64) = params(63);

%% Params
try
    params(87:101)     = [0,0,1,0,0,0,0,0,1,0,0,0,1,0,0];  % Steady State Meals are solid and there is no meal in t0
    
    params(102:105)   = [0 , 1, 180, 75.6];                % female, male, height, weight
    params(106:109)   = [1,0,0,0];                         % This is meal A
    params(110:117)   = [10, 0, 138.64, 29.25,1,0,0,0];    % Meal Information
    simTaylorA  = model(time,[], params);
    
    params(106:109)   = [0,1,0,0];                         % This is meal B
    simTaylorB  = model(time,[], params);
    
catch error
    
    Cost_Model = 1e60;
    return
end


%% Summerize Cost
simG          = simTaylorB.variablevalues(:,ismember(simTaylorB.variables,'G'));
simI          = simTaylorB.variablevalues(:,ismember(simTaylorB.variables,'I'));
simGly        = simTaylorA.variablevalues(:,ismember(simTaylorA.variables,'Glycogen_liver'));
    
simGly = simGly(round(Taylor1996_data.time_glycogen_healthy(1:9),0));
simG   = simG(round(Taylor1996_data.time_glucose_healthy(1:13),0));
simI   = simI(round(Taylor1996_data.time_Insulin_healthy(1:13),0));

CostGly     = (Taylor1996_data.glycogen_healthy(1:datapointsGly) - simGly(1:datapointsGly)).^2./(Taylor1996_data.glycogenSEM_healthy(1:datapointsGly).^2);
CostG       = ((Taylor1996_data.glucose_healthy(1:datapointsG) - simG(1:datapointsG)).^2)  ./(Taylor1996_data.glucoseSEM_healthy(1:datapointsG).^2);
CostI       = ((Taylor1996_data.Insulin_healthy(1:datapointsI) - simI(1:datapointsI)).^2)  ./(Taylor1996_data.InsulinSEM_healthy(1:datapointsI).^2);

Cost_Model = nansum(CostGly,"all") + nansum(CostG,"all") + nansum(CostI,"all");

end

