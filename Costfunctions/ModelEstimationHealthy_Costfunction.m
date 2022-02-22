function Cost_Model = ModelEstimationHealthy_Costfunction(Krssak2004,Lerche2009,Magnusson1992,Silfvergren2021,time,params,modelName)

model = str2func(modelName);
params=exp(params);

%% Try Simulation

try
    params(99:101)    = [1,0,0];                         % Steady State Meals are solid and there is no meal in t0
    
    params(102:105)   = [0 , 1, 175, 75];                % female, male, height, weight
    params(106:109)   = [1,0,0,0];                       % This is study A
    params(110:117)   = [10, 0, 87, 23,1,0,0,0];         % Meal Information
    simKrssak2004  = model(time,[], params);
    
    params(102:105)   = [0 , 1, 175, 75];                % female, male, height, weight
    params(106:109)   = [0,1,0,0];                       % This is study B
    params(110:117)   = [10, 2420, 75, 0, 1, 0,0,0];     % Meal Information
    simLerche2009  = model(time,[], params);
    
    params(102:105)    = [0 , 1, 175, 75 ];               % female, male, height, weight
    params(106:109)    = [0,0,1,0];                       % This is study C
    params(110:117)    = [5, 0, 98.3125, 26, 1, 0,0,0];   % Meal Information
    simMagnusson1992   = model(time,[], params);
    
    params(102:105)     = [0 , 1, 175, 75 ];                    % female, male, height, weight
    params(106:109)     = [0,0,0,1];                            % This is study D
    params(110:117)     = [10, 0, 155.35, 20.7, 1, 0,0,0];   % Meal Information
    simSilfvergren2021  = model(time,[], params);
    
catch error
    simKrssak2004      = ones(size(time))*10^5;
    simMagnusson1992   = ones(size(time))*10^5;
    simLerche2009      = ones(size(time))*10^5;
    simSilfvergren2021 = ones(size(time))*10^5;
    
    Cost_Model = 1e60;
    return
end

%% Krssak
simGly = simKrssak2004.variablevalues(:,ismember(simKrssak2004.variables,'Glycogen_liver'));
simG   = simKrssak2004.variablevalues(:,ismember(simKrssak2004.variables,'G'));
simI   = simKrssak2004.variablevalues(:,ismember(simKrssak2004.variables,'I'));
simEGP = simKrssak2004.reactionvalues(:,ismember(simKrssak2004.reactions,'EGP'));

simGly = simGly(round(Krssak2004.timeGlycogen(1:22),0));
simG   = simG(round(Krssak2004.timeGlucose(1:25),0));
simI   = simI(round(Krssak2004.timeInsulin(1:10),0));
simEGP = simEGP(round(Krssak2004.timeEGP(1:20),0));

CostEGP = ((Krssak2004.EGP(1:20)     - simEGP).^2)./((Krssak2004.EGPSEM(1:20)).^2);
CostGly = ((Krssak2004.glycogen(1:22)- simGly).^2)./(Krssak2004.glycogenSEM(1:22).^2);
CostG   = ((Krssak2004.glucose(1:25) - simG).^2)  ./(Krssak2004.glucoseSEM(1:25).^2);
CostI   = ((Krssak2004.insulin(1:10) - simI).^2)  ./((Krssak2004.insulinSEM(1:10)).^2);


Cost_Model = nansum(CostG,"all") + nansum(CostI,"all") + nansum(CostGly,"all") + nansum(CostEGP,"all");

%% Silfvergren 2021

simG   = simSilfvergren2021.variablevalues(:,ismember(simSilfvergren2021.variables,'G'))/18;
simG   = simG(round(Silfvergren2021.time(1:30),0));
CostG  = ((Silfvergren2021.value(1:30) - simG).^2)  ./((Silfvergren2021.value(1:30)*0.15).^2);

Cost_Model = Cost_Model + nansum(CostG,"all");

%% Lerche2009
simG   = simLerche2009.variablevalues(:,ismember(simLerche2009.variables,'G'));
simI   = simLerche2009.variablevalues(:,ismember(simLerche2009.variables,'I'));
simG   = simG(round(Lerche2009.time_glucose_healthy(1:56),0));
simI   = simI(round(Lerche2009.time_insulin_healthy(1:21),0));

CostG = (Lerche2009.glucose_healthy(1:56) - simG).^2       ./(Lerche2009.glucoseSEM_healthy(1:56).^2);
CostI = (Lerche2009.insulin_healthy(1:21) - simI(1:21)).^2 ./(Lerche2009.insulinSEMFIXED_healthy(1:21).^2);

Cost_Model = Cost_Model + nansum(CostG,"all") + nansum(CostI,"all");

%% Magnusson1992

simGly = simMagnusson1992.variablevalues(:,ismember(simMagnusson1992.variables,'Glycogen_liver'));
simGly = simGly(round(Magnusson1992.time_glycogen_healthy(1:5),0));
CostGly = (Magnusson1992.glycogen_healthy(1:5) - simGly).^2./((Magnusson1992.glycogenSEM_healthy(1:5)).^2);

Gluconeogenesis    = sum(simMagnusson1992.reactionvalues(7914:9294,ismember(simMagnusson1992.reactions,'Gluconeogenesis'))) + sum(simMagnusson1992.reactionvalues(7914:9294,ismember(simMagnusson1992.reactions,'EGP_Kidneys')));
EGP                = sum(simMagnusson1992.reactionvalues(7914:9294,ismember(simMagnusson1992.reactions,'EGP')));
Procentage         = Gluconeogenesis/EGP;
CostEGP = ((0.70 - Procentage)^2)/(0.06^2);

Cost_Model = Cost_Model + nansum(CostGly,"all") + nansum(CostEGP,"all");

if Procentage > 1
    Cost_Model = Cost_Model + 500;
end

end
