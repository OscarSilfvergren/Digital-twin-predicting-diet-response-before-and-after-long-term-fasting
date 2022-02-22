function Cost_Model = ModelEstimationT2DM_Costfunction(Krssak2004_Diabetes,Magnusson1992,time,params,modelName)

model = str2func(modelName);
params=exp(params);

%% Try Simulation

try
    params(99:101)  = [1,0,0];                           % Steady State Meals are solid and there is no meal in t0
    
    params(102:105)   = [0 , 1, 175, 75];                % female, male, height, weight
    params(106:109)   = [1,0,0,0];                       % This is study A
    params(110:117)   = [10, 0, 87, 23,1,0,0,0];         % Meal Information
    simKrssak2004  = model(time,[], params);
    
    params(102:105)    = [0 , 1, 175, 75 ];               % female, male, height, weight
    params(106:109)    = [0,0,1,0];                       % This is study C
    params(110:117)    = [5, 0, 98.3125, 26, 1, 0,0,0];   % Meal Information
    simMagnusson1992   = model(time,[], params);
    
catch error
    simKrssak2004     = ones(size(time))*10^5;
    simMagnusson1992  = ones(size(time))*10^5;
    Cost_Model = 1e60;
    return
end

%% Krssak 2004

simGly = simKrssak2004.variablevalues(:,ismember(simKrssak2004.variables,'Glycogen_liver'));
simG   = simKrssak2004.variablevalues(:,ismember(simKrssak2004.variables,'G'));
simI   = simKrssak2004.variablevalues(:,ismember(simKrssak2004.variables,'I'));
simEGP = simKrssak2004.reactionvalues(:,ismember(simKrssak2004.reactions,'EGP'));

simGly = simGly(round(Krssak2004_Diabetes.timeGlycogen(1:22),0));
simG   = simG(round(Krssak2004_Diabetes.timeGlucose(1:26),0));
simI   = simI(round(Krssak2004_Diabetes.timeInsulin(1:10),0));
simEGP = simEGP(round(Krssak2004_Diabetes.timeEGP(1:20),0));

CostGly = ((Krssak2004_Diabetes.glycogen(1:22)- simGly).^2)./(Krssak2004_Diabetes.glycogenSEM(1:22).^2);
CostG   = ((Krssak2004_Diabetes.glucose(1:26) - simG).^2)  ./(Krssak2004_Diabetes.glucoseSEM(1:26).^2);
CostI   = ((Krssak2004_Diabetes.insulin(1:10) - simI).^2)  ./((Krssak2004_Diabetes.insulinSEM(1:10)).^2);
CostEGP = ((Krssak2004_Diabetes.EGP(1:20)     - simEGP).^2)./(Krssak2004_Diabetes.EGPSEM(1:20).^2);

Cost_Model =  nansum(CostG,"all") + nansum(CostI,"all") + nansum(CostGly,"all") + nansum(CostEGP,"all");

%% Magnusson 1992

simGly  = simMagnusson1992.variablevalues(:,ismember(simMagnusson1992.variables,'Glycogen_liver'));
simGly = simGly(round(Magnusson1992.time_glycogen_diabetes(1:5),0));
CostGly = (Magnusson1992.glycogen_diabetes(1:5) - simGly).^2./((Magnusson1992.glycogenSEM_diabetes(1:5)).^2);

Gluconeogenesis    = sum(simMagnusson1992.reactionvalues(7914:9294,ismember(simMagnusson1992.reactions,'Gluconeogenesis'))) + sum(simMagnusson1992.reactionvalues(7914:9294,ismember(simMagnusson1992.reactions,'EGP_Kidneys')));
EGP                = sum(simMagnusson1992.reactionvalues(7914:9294,ismember(simMagnusson1992.reactions,'EGP')));
Procentage         = Gluconeogenesis/EGP;
CostEGP = ((0.88 - Procentage)^2)/(0.02^2);

Cost_Model = Cost_Model + nansum(CostGly,"all") + nansum(CostEGP,"all");
end
