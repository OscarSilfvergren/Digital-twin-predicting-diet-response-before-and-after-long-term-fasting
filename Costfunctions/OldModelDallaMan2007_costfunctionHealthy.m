function Cost_Model = OldModelDallaMan2007_costfunctionHealthy(DallaMan2007_data,params,modelName)

model = str2func(modelName);
params=exp(params);

params(41:43) = [0,60,75000];
time = linspace(0,8500,8501);


%% Try Simulation

try
    sim = model(time,[], params);
catch error
    sim = ones(size(time))*10^5;
    Cost_Model = 1e60;
    return
end

%% Fitting to data
simG    = sim.reactionvalues(:,ismember(sim.reactions,'plasmaGlucose'));
simI    = sim.reactionvalues(:,ismember(sim.reactions,'plasmaInsulin'));
simEGP  = sim.reactionvalues(:,ismember(sim.reactions,'endogenousGlucoseProduction'));
simRa   = sim.reactionvalues(:,ismember(sim.reactions,'glucoseRateOfAppearance'));
simU    = sim.reactionvalues(:,ismember(sim.reactions,'glucoseUptake'));
simS    = sim.reactionvalues(:,ismember(sim.reactions,'insulinSecretion'));

simG   = simG(round(DallaMan2007_data.Glucose_time(1:21),0)+7680);
simI   = simI(round(DallaMan2007_data.Insulin_time(1:21),0)+7680);
simEGP = simEGP(round(DallaMan2007_data.EGP_time(1:22),0)+7680);
simRa  = simRa(round(DallaMan2007_data.Ra_time(1:20),0)+7680);
simU   = simU(round(DallaMan2007_data.U_time(1:21),0)+7680);
simS   = simS(round(DallaMan2007_data.S_time(1:21),0)+7680);

CostG   = (DallaMan2007_data.Glucose(1:21) - simG).^2./(DallaMan2007_data.GlucoseSEM(1:21).^2);
CostI   = (DallaMan2007_data.Insulin(1:21) - simI).^2./(DallaMan2007_data.InsulinSEM(1:21).^2);
CostEGP = (DallaMan2007_data.EGP(1:22) - simEGP).^2./(DallaMan2007_data.EGPSEM(1:22).^2);
CostRa  = (DallaMan2007_data.Ra(1:20) - simRa).^2./(DallaMan2007_data.RaSEM(1:20).^2);
CostU   = (DallaMan2007_data.U(1:21) - simU).^2./(DallaMan2007_data.USEM(1:21).^2);
CostS   = (DallaMan2007_data.S(1:21) - simS).^2./(DallaMan2007_data.SSEM(1:21).^2);

Cost_Model =  nansum(CostRa,"all") + nansum(CostG,"all") + nansum(CostI,"all") + nansum(CostEGP,"all") + nansum(CostU,"all") + nansum(CostS,"all");

end

