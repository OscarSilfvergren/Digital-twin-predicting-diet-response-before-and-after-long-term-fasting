function Cost_Model = Silfvergren2021_costfunction(data,time_data,FedState,time,params,modelName,body_information,meal_information,datapointsG)

model = str2func(modelName);
params=exp(params);

params  = AssignParameter(params, body_information, meal_information);

% Simulation
try
    sim = model(time,[],params);
    
catch error
    sim = ones(size(time))*10^5;
    Cost_Model = 1e60;
    return
end


%% Create weight to data without SEM

simG   = sim.variablevalues(:,ismember(sim.variables,'G'))/18;

simulation   = simG(time_data);

GlucoseWeight    = data(1:datapointsG)*0.1;
Cost   = ((data(1:datapointsG) - simulation(1:datapointsG)).^2)  ./((GlucoseWeight).^2); %Note that "GlucoseWeight" is used as no SEM is aviable

Cost_Model = nansum(Cost,"all");

%% Temp
if FedState ==1
    simGly =   sim.variablevalues(time_data(1),ismember(sim.variables,'Glycogen_liver'));
    
    if simGly < 250
    Cost_Model = Cost_Model + 1000;
    end
    
elseif FedState == 0
    simGly =   sim.variablevalues(time_data(1),ismember(sim.variables,'Glycogen_liver'));
    
        if simGly > 50
    Cost_Model = Cost_Model + 1000;
    return
        end
    
end

end

