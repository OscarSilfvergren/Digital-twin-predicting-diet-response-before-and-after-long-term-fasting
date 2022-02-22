function Cost_Model = Silfvergren2021_costfunction(data,time_data,time,params,modelName,body_information,meal_information,datapointsG)

model = str2func(modelName);
params=exp(params);

params  = AssignParameter(params, body_information, meal_information);

try
    sim = model(time,[],params);
    
catch error
    sim = ones(size(time))*10^5;
    Cost_Model = 1e60;
    return
end


%% p1
simG   = sim.variablevalues(:,ismember(sim.variables,'G'))/18;

simulation   = simG(time_data);

SEM    = data(1:datapointsG)*0.1;
Cost   = ((data(1:datapointsG) - simulation(1:datapointsG)).^2)  ./((SEM).^2);

Cost_Model = nansum(Cost,"all");

end

