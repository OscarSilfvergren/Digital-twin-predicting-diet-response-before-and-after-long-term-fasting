function Cost_Model = Firth1986_costfunctionHealthy(Realvalues,time,params,modelName,meal_information,body_information,datapointEGP,datapointG,datapointI);

model = str2func(modelName);
params=exp(params);

EGP         = Realvalues.EGP_healthy(1:14);
glucose     = Realvalues.glucose_healthy(1:15);
insulin     = Realvalues.insulin_healthy(1:15);

EGP_SEM     = Realvalues.EGPSEM_healthy(1:14);
glucose_SEM = Realvalues.glucoseSEM_healthy(1:15);
insulin_SEM = Realvalues.insulinSEM_healthy(1:15);

params = AssignParameter(params, body_information, meal_information);

try
    sim = model(time,[], params);
catch error
    sim = ones(size(time))*10^5;
    Cost_Model = 1e60;
    return
end


%% Summerize Cost
simEGP = sim.reactionvalues(:,ismember(sim.reactions,'EGP'));
simG   = sim.variablevalues(:,ismember(sim.variables,'G'));
simI   = sim.variablevalues(:,ismember(sim.variables,'I'));

simG   = simG(round(Realvalues.time_glucose_healthy(1:15),0));
simI   = simI(round(Realvalues.time_insulin_healthy(1:15),0));
simEGP = simEGP(round(Realvalues.time_EGP_healthy(1:14),0));

CostEGP = (EGP(1:datapointEGP) - simEGP(1:datapointEGP)).^2./(EGP_SEM(1:datapointEGP).^2);
CostG   = (glucose(1:datapointG) - simG(1:datapointG)).^2./(glucose_SEM(1:datapointG).^2);
CostI   = (insulin(1:datapointI) - simI(1:datapointI)).^2./(insulin_SEM(1:datapointI).^2);

Cost_Model =  nansum(CostG,"all") + nansum(CostI,"all") + nansum(CostEGP,"all");

end
