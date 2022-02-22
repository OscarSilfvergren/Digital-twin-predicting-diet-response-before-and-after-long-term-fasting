function Cost_Model = Rothman1991_costfunctionIndividual(data,data_time,time,params,modelName,body_information, meal_information)

model = str2func(modelName);
params=exp(params);

params = AssignParameter(params, body_information, meal_information);

try
    sim = model(time,[], params);
catch error
    sim = ones(size(time))*10^5;
    Cost_Model = 1e60;
    return
end


%% Data
GlycogenSEM = (data*0.1)+20;

simGly      = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
simGly      = simGly(round(data_time,0));
CostGly     = (data - simGly).^2./(GlycogenSEM.^2);

Cost_Model = nansum(CostGly,"all");

end

