function Cost_Model = EstimationData_costfunction(Krssak2004Healthy_data,Rothman1991_data,Silfvergren2021_data,params,model,initialvalues, pNames,sNames,AmountParametersOptimized)

params=exp(params);

%% Declare individual params

% Study 2
Study2_SteadystateMeal         = [params(AmountParametersOptimized+1) , params(AmountParametersOptimized+2)]; % Carbs and Protein flow
Study2_BasalInsulin            = params(AmountParametersOptimized+3); % BasalInsulin
Study2_BasalGlucose            = params(AmountParametersOptimized+4); % BasalGlucose
Study2_InsulinProduction       = params(AmountParametersOptimized+5); % InsulinProduction
Study2_InsulinClearance        = params(AmountParametersOptimized+6); % InsulinClearance
Study2_InsulinResponse         = params(AmountParametersOptimized+7); % InsulinResponse
Study2_BloodVolumeUncertainty  = params(AmountParametersOptimized+8); % BloodVolumeUncertainty
Study2_BloodLiverUncerteinty   = params(AmountParametersOptimized+9); % BloodLiverUncertainty

% Study 3
Study3_SteadystateMeal          = [params(AmountParametersOptimized+10) , params(AmountParametersOptimized+11)]; % Carbs and Protein flow
Study3_BasalInsulin             = params(AmountParametersOptimized+12); % BasalInsulin
Study3_BasalGlucose             = params(AmountParametersOptimized+13); % BasalGlucose
Study3_InsulinProduction        = params(AmountParametersOptimized+14); % InsulinProduction
Study3_InsulinClearance         = params(AmountParametersOptimized+15); % InsulinClearance
Study3_InsulinResponse          = params(AmountParametersOptimized+16); % InsulinResponse
Study3_BloodVolumeUncertainty   = params(AmountParametersOptimized+17);  % BloodVolumeUncertainty
Study3_BloodLiverUncerteinty    = params(AmountParametersOptimized+18);  % BloodLiverUncertainty


%% Simulate study 1 - Krssak2004

params(ismember(pNames,'meal_solid'))     = 1;
params(ismember(pNames,'meal_liquid'))    = 0;
params(ismember(pNames,'D'))              = params(ismember(pNames,'Meal_CarbohydratesFlow'))*1440; % mg
params(ismember(pNames,'female_boolean')) = 0;
params(ismember(pNames,'male_boolean'))   = 1;
params(ismember(pNames,'Height'))         = 175; % cm
params(ismember(pNames,'BW_start'))       = 75;  % kg

paramsKrsakkMeal = params;
paramsKrsakkMeal(ismember(pNames,'Meal_CarbohydratesFlow')) = 87*1000/15;   % mg/min
paramsKrsakkMeal(ismember(pNames,'Meal_ProteinFlow'))       = 23*1000/15;   % mg/min
paramsKrsakkMeal(ismember(pNames,'D'))                      = 87*1000;      % mg

try
    
    SimSteadystate       = model(0:40320:40320,initialvalues,params(1:65)); % steady state sim with constant food (1 month)
    params(ismember(pNames,"Meal_CarbohydratesFlow"))    = 0;
    params(ismember(pNames,"Meal_ProteinFlow"))          = 0;
    SimSteadystate       = model(0:360:360,SimSteadystate.statevalues(end,:),params(1:65)); % steady state sim with no food
    
    SimKrssakToMeal                = model(1:1:360,SimSteadystate.statevalues(end,:),params(1:65));                % Sim to meal
    SimKrssak_meal                 = model(360:1:375,SimKrssakToMeal.statevalues(end,:),paramsKrsakkMeal(1:65));   % Sim  meal consumtion
    params(ismember(pNames,'D'))   = paramsKrsakkMeal(ismember(pNames,'D'));                                       % emtying of stomach is based on new meal not steady state meals
    SimKrssakAfterMeal             = model(375:1:1260,SimKrssak_meal.statevalues(end,:),params(1:65));             % Sim  after meal
    
catch error
    Cost_Model = 1e60;
    return
end


% Merge
simKrssak2004.time           = 1:1:1260;
simKrssak2004.variablevalues = [SimKrssakToMeal.variablevalues(1:end-1,:); SimKrssak_meal.variablevalues(1:end-1,:); SimKrssakAfterMeal.variablevalues(1:end,:)];
simKrssak2004.reactionvalues = [SimKrssakToMeal.reactionvalues(1:end-1,:); SimKrssak_meal.reactionvalues(1:end-1,:); SimKrssakAfterMeal.reactionvalues(1:end,:)];
simKrssak2004.statevalues    = [SimKrssakToMeal.statevalues(1:end-1,:);    SimKrssak_meal.statevalues(1:end-1,:);    SimKrssakAfterMeal.statevalues(1:end,:)];
simKrssak2004.variables      = [SimKrssakToMeal.variables];
simKrssak2004.reactions      = [SimKrssakToMeal.reactions];
simKrssak2004.states         = [SimKrssakToMeal.states];

simGly = simKrssak2004.variablevalues(:,ismember(simKrssak2004.variables,'Glycogen_liver'));
simG   = simKrssak2004.variablevalues(:,ismember(simKrssak2004.variables,'G'));
simI   = simKrssak2004.variablevalues(:,ismember(simKrssak2004.variables,'I'));
simEGP = simKrssak2004.reactionvalues(:,ismember(simKrssak2004.reactions,'EGP'));

simGly = simGly(round(Krssak2004Healthy_data.timeGlycogen(1:22)+360,0)); % +360 because meal is at 360 min into simulation
simG   = simG(round(Krssak2004Healthy_data.timeGlucose(1:25)+360,0));
simI   = simI(round(Krssak2004Healthy_data.timeInsulin(1:10)+360,0));
simEGP = simEGP(round(Krssak2004Healthy_data.timeEGP(1:20)+360,0));

CostEGP = ((Krssak2004Healthy_data.EGP(1:20)     - simEGP).^2)./((Krssak2004Healthy_data.EGPSEM(1:20)).^2);
CostGly = ((Krssak2004Healthy_data.glycogen(1:22)- simGly).^2)./(Krssak2004Healthy_data.glycogenSEM(1:22).^2);
CostG   = ((Krssak2004Healthy_data.glucose(1:25) - simG).^2)  ./(Krssak2004Healthy_data.glucoseSEM(1:25).^2);
CostI   = ((Krssak2004Healthy_data.insulin(1:10) - simI).^2)  ./((Krssak2004Healthy_data.insulinSEM(1:10)).^2);

Cost_Model = nansum(CostG,"all") + nansum(CostI,"all") + nansum(CostGly,"all") + nansum(CostEGP,"all");

%% Qualitative demands to ensure reasonable organ specific behaviour

% fed state organ specific glucose produciton
EGP_Kidneys = sum(simKrssak2004.reactionvalues(360:720,ismember(simKrssak2004.reactions,'EGP_Kidneys')));
EGP_Liver   = sum(simKrssak2004.reactionvalues(360:720,ismember(simKrssak2004.reactions,'EGP_Liver')));
EGP         = sum(simKrssak2004.reactionvalues(360:720,ismember(simKrssak2004.reactions,'EGP')));
ratioK       = EGP_Kidneys/EGP; % circa 20% +- 15%
ratioL       = EGP_Liver/EGP;   % circa 80% +- 15%

%ratioK
if ratioK > 0.35
    Cost_Model = Cost_Model + 5000;
end
if 0.05 > ratioK
    Cost_Model = Cost_Model + 5000;
end

%ratioL
if ratioL > 0.95
    Cost_Model = Cost_Model + 5000;
end
if 0.45 > ratioL
    Cost_Model = Cost_Model + 5000;
end

% Glucose specific organ utilization 6h before meal
U_idl = sum(SimKrssakToMeal.reactionvalues(:,ismember(SimKrssakToMeal.reactions,'U_il')));
U_idm = sum(SimKrssakToMeal.reactionvalues(:,ismember(SimKrssakToMeal.reactions,'U_idm')));
U_idf = sum(SimKrssakToMeal.reactionvalues(:,ismember(SimKrssakToMeal.reactions,'U_idf')));
U_ii  = sum(SimKrssakToMeal.reactionvalues(:,ismember(SimKrssakToMeal.reactions,'U_ii_print')));
U     = sum(SimKrssakToMeal.reactionvalues(:,ismember(SimKrssakToMeal.reactions,'U')));
ratio_f = U_idf/U; % circa 17.5% +- 10%
ratio_m = U_idm/U; % circa 22.5% +- 10%
ratio_b = U_ii/U;  % circa 40% +- 10%
ratio_l = U_idl/U; % circa 20% +- 10%

%ratio_f
if ratio_f > 0.275
    Cost_Model = Cost_Model + 5000;
end
if 0.075 > ratio_f
    Cost_Model = Cost_Model + 5000;
end

%ratio_m
if ratio_m > 0.325
    Cost_Model = Cost_Model + 5000;
end
if 0.125 > ratio_m
    Cost_Model = Cost_Model + 5000;
end

%ratio_b
if ratio_b > 0.50
    Cost_Model = Cost_Model + 5000;
end
if 0.30 > ratio_b
    Cost_Model = Cost_Model + 5000;
end

%ratio_l
if ratio_l > 0.30
    Cost_Model = Cost_Model + 5000;
end
if 0.10 > ratio_l
    Cost_Model = Cost_Model + 5000;
end

% Insulin clearance 6h starting at meal
InsulinDegradation_Blood     = sum(simKrssak2004.reactionvalues(360:720,ismember(simKrssak2004.reactions,'InsulinDegradation_Blood')));
InsulinDegradation_Liver     = sum(simKrssak2004.reactionvalues(360:720,ismember(simKrssak2004.reactions,'InsulinDegradation_Liver')));
ratioLiver = InsulinDegradation_Liver/(InsulinDegradation_Blood + InsulinDegradation_Liver); % circa 75% +- 10%
ratioBlood = InsulinDegradation_Blood/(InsulinDegradation_Blood + InsulinDegradation_Liver); % circa 25% +- 10%

%ratioLiver
if ratioLiver > 0.85
    Cost_Model = Cost_Model + 5000;
end
if 0.65 > ratioLiver
    Cost_Model = Cost_Model + 5000;
end

%ratioBlood
if ratioBlood > 0.35
    Cost_Model = Cost_Model + 5000;
end
if 0.15 > ratioBlood
    Cost_Model = Cost_Model + 5000;
end

%% Simulate study 2 - Rothman1991_data

% Declare personal scaling
%Meals
params(ismember(pNames,"Meal_CarbohydratesFlow"))    = Study2_SteadystateMeal(1); % Carbs
params(ismember(pNames,"Meal_ProteinFlow"))          = Study2_SteadystateMeal(2); % Protein
% Basal Values: Personal Scaling -  Param 3 & 4
params(ismember(pNames,"S_b"))          = params(ismember(pNames,"S_b")) * Study2_BasalInsulin;
params(ismember(pNames,"G_b"))          = params(ismember(pNames,"S_b")) * Study2_BasalGlucose;
% InsulinProduction: Personal Scaling - Param 5
params(ismember(pNames,"alpha"))          = params(ismember(pNames,"alpha")) * Study2_InsulinProduction;
params(ismember(pNames,"beta"))           = params(ismember(pNames,"beta"))  * Study2_InsulinProduction;
% InsulinClearance: Personal Scaling -  Param 6
params(ismember(pNames,"m_3"))          = params(ismember(pNames,"m_3"))   *Study2_InsulinClearance;
params(ismember(pNames,"m_4"))          = params(ismember(pNames,"m_4"))   *Study2_InsulinClearance;
% InsulinResponse: Personal Scaling -  Param 7
params(ismember(pNames,"V_mX"))          = params(ismember(pNames,"V_mX"))  * Study2_InsulinResponse;
params(ismember(pNames,"V_fX"))          = params(ismember(pNames,"V_fX"))  * Study2_InsulinResponse;
% BloodVolumeUncertainty: Personal Scaling -  Param 8 & 9
params(ismember(pNames,"BloodVolumeUncertainty"))          = params(ismember(pNames,"BloodVolumeUncertainty"))  * Study2_BloodVolumeUncertainty;
params(ismember(pNames,"BloodLiverUncerteinty"))           = params(ismember(pNames,"BloodLiverUncerteinty"))   * Study2_BloodLiverUncerteinty;

params(ismember(pNames,'meal_solid'))     = 1;
params(ismember(pNames,'meal_liquid'))    = 0;
params(ismember(pNames,'D'))              = params(ismember(pNames,'Meal_CarbohydratesFlow'))*1440; % mg
params(ismember(pNames,'female_boolean')) = 0;
params(ismember(pNames,'male_boolean'))   = 1;
params(ismember(pNames,'Height'))         = 180; %cm
params(ismember(pNames,'BW_start'))       = 75;  %kg

paramsRothman = params;
paramsRothman(ismember(pNames,'meal_CarbohydratesFlow')) = 0;      % mg/min - no meals in study
paramsRothman(ismember(pNames,'meal_ProteinFlow'))       = 0;      % mg/min
paramsRothman(ismember(pNames,'D'))                      = 0;      % mg

try
    SimSteadystate       = model(0:40320:40320,initialvalues,params(1:65));            % steady state sim meals
    params(ismember(pNames,"Meal_CarbohydratesFlow"))  = 0;
    params(ismember(pNames,"Meal_ProteinFlow"))        = 0;
    SimSteadystate       = model(0:360:360,SimSteadystate.statevalues(end,:),params(1:65)); % steady state sim no meals
    
    SimRothman      = model(1:1:3700,SimSteadystate.statevalues(end,:),paramsRothman(1:65)); % study sim based on steady state sim
    
    
catch error
    Cost_Model = 1e60;
    return
end

simGly  = SimRothman.variablevalues(:,ismember(SimRothman.variables,'Glycogen_liver'));
simGly  = simGly(round(Rothman1991_data.timepALL(1:11),0));
CostGly = ((Rothman1991_data.glycogenpALL(1:11) - simGly(1:11)).^2)./(Rothman1991_data.glycogenpSEM(1:11).^2);

Gluconeogenesis1    = sum(SimRothman.reactionvalues(1:1320,ismember(SimRothman.reactions,'Gluconeogenesis')));
Gluconeogenesis2    = sum(SimRothman.reactionvalues(1320:2160,ismember(SimRothman.reactions,'Gluconeogenesis')));
Gluconeogenesis3    = sum(SimRothman.reactionvalues(2160:3240,ismember(SimRothman.reactions,'Gluconeogenesis')));
EGP1                = sum(SimRothman.reactionvalues(1:1320,ismember(SimRothman.reactions,'EGP')));
EGP2                = sum(SimRothman.reactionvalues(1320:2160,ismember(SimRothman.reactions,'EGP')));
EGP3                = sum(SimRothman.reactionvalues(2160:3240,ismember(SimRothman.reactions,'EGP')));
Procentage1         = Gluconeogenesis1/EGP1;
Procentage2         = Gluconeogenesis2/EGP2;
Procentage3         = Gluconeogenesis3/EGP3;

% GNG can not be more than 100% of EGP
if Procentage3 > 1
    Cost_Model = Cost_Model + 5000;
    
end

CostGNG = (([0.64, 0.82,  0.96] - [Procentage1, Procentage2, Procentage3]).^2)./([0.05, 0.05, 0.01].^2);

Cost_Model = Cost_Model + nansum(CostGNG,"all") + nansum(CostGly,"all");

%% Simulate study 3 - Silfvergren2021_data

% Declare personal scaling
%Meals
params(ismember(pNames,"Meal_CarbohydratesFlow"))    = Study3_SteadystateMeal(1); % Carbs
params(ismember(pNames,"Meal_ProteinFlow"))          = Study3_SteadystateMeal(2); % Protein
% Basal Values: Personal Scaling -  Param 3 & 4
params(ismember(pNames,"S_b"))          = params(ismember(pNames,"S_b")) * Study3_BasalInsulin;
params(ismember(pNames,"G_b"))          = params(ismember(pNames,"S_b")) * Study3_BasalGlucose;
% InsulinProduction: Personal Scaling - Param 5
params(ismember(pNames,"alpha"))          = params(ismember(pNames,"alpha")) * Study3_InsulinProduction;
params(ismember(pNames,"beta"))           = params(ismember(pNames,"beta"))  * Study3_InsulinProduction;
% InsulinClearance: Personal Scaling -  Param 6
params(ismember(pNames,"m_3"))          = params(ismember(pNames,"m_3"))   *Study3_InsulinClearance;
params(ismember(pNames,"m_4"))          = params(ismember(pNames,"m_4"))   *Study3_InsulinClearance;
% InsulinResponse: Personal Scaling -  Param 7
params(ismember(pNames,"V_mX"))          = params(ismember(pNames,"V_mX"))  * Study3_InsulinResponse;
params(ismember(pNames,"V_fX"))          = params(ismember(pNames,"V_fX"))  * Study3_InsulinResponse;
% BloodVolumeUncertainty: Personal Scaling -  Param 8 & 9
params(ismember(pNames,"BloodVolumeUncertainty"))          = params(ismember(pNames,"BloodVolumeUncertainty"))  * Study3_BloodVolumeUncertainty;
params(ismember(pNames,"BloodLiverUncerteinty"))           = params(ismember(pNames,"BloodLiverUncerteinty"))   * Study3_BloodLiverUncerteinty;

params(ismember(pNames,'meal_solid'))     = 1;
params(ismember(pNames,'meal_liquid'))    = 0;
params(ismember(pNames,'D'))              = params(ismember(pNames,"Meal_CarbohydratesFlow"))*1440; % mg
params(ismember(pNames,'female_boolean')) = 0;
params(ismember(pNames,'male_boolean'))   = 1;
params(ismember(pNames,'Height'))         = 178; %cm
params(ismember(pNames,'BW_start'))       = 80;  %kg

paramsSilfvergrenMeal                                            = params;
paramsSilfvergrenMeal(ismember(pNames,'Meal_ProteinFlow'))       = 25.55*1000/3;  % mg/min
paramsSilfvergrenMeal(ismember(pNames,'Meal_CarbohydratesFlow')) = 2.6*1000/3;    % mg/min
paramsSilfvergrenMeal(ismember(pNames,'D'))                      = 2.6*1000;      % mg
paramsSilfvergrenMeal(ismember(pNames,'meal_solid'))     = 1;
paramsSilfvergrenMeal(ismember(pNames,'meal_liquid'))    = 0;

try
    
    SimSteadystate       = model(0:40320:40320,initialvalues,params(1:65));
    params(ismember(pNames,"Meal_CarbohydratesFlow"))    = 0;
    params(ismember(pNames,"Meal_ProteinFlow"))          = 0;
    SimSteadystate       = model(0:360:360,SimSteadystate.statevalues(end,:),params(1:65));
    
    SimSilfvergrenToFedMeal       = model(1:1:19,SimSteadystate.statevalues(end,:),params(1:65));                           % short sim to fed state meal
    SimSilfvergrenFedMeal         = model(19:1:22,SimSilfvergrenToFedMeal.statevalues(end,:),paramsSilfvergrenMeal(1:65));  % consumtion of fed state meal
    params(ismember(pNames,'D'))  = 2.6*1000;   % mg
    SimSilfvergrenToUnfedMeal     = model(22:1:2881,SimSilfvergrenFedMeal.statevalues(end,:),params(1:65));                 % sim to unfed meal
    SimSilfvergrenUnfedMeal       = model(2881:1:2884,SimSilfvergrenToUnfedMeal.statevalues(end,:),paramsSilfvergrenMeal(1:65));  % consumtion of unfed state meal
    SimSilfvergrenAfterUnfedMeal  = model(2884:1:3300,SimSilfvergrenUnfedMeal.statevalues(end,:),params(1:65));                   % sim after unfed meal
    
catch error
    Cost_Model = 1e60;
    return
end

% Merge
simSilfvergren.time           = 1:1:3300;
simSilfvergren.variablevalues = [SimSilfvergrenToFedMeal.variablevalues(1:end-1,:); SimSilfvergrenFedMeal.variablevalues(1:end-1,:); SimSilfvergrenToUnfedMeal.variablevalues(1:end-1,:); SimSilfvergrenUnfedMeal.variablevalues(1:end-1,:); SimSilfvergrenAfterUnfedMeal.variablevalues(1:end,:)];
simSilfvergren.reactionvalues = [SimSilfvergrenToFedMeal.reactionvalues(1:end-1,:); SimSilfvergrenFedMeal.reactionvalues(1:end-1,:); SimSilfvergrenToUnfedMeal.reactionvalues(1:end-1,:); SimSilfvergrenUnfedMeal.reactionvalues(1:end-1,:); SimSilfvergrenAfterUnfedMeal.reactionvalues(1:end,:)];
simSilfvergren.statevalues    = [SimSilfvergrenToFedMeal.statevalues              ; SimSilfvergrenFedMeal.statevalues(1:end-1,:)   ; SimSilfvergrenToUnfedMeal.statevalues(1:end-1,:)   ; SimSilfvergrenUnfedMeal.statevalues(1:end-1,:)   ; SimSilfvergrenAfterUnfedMeal.statevalues(1:end,:)];
simSilfvergren.variables      = [SimSilfvergrenToFedMeal.variables];
simSilfvergren.reactions      = [SimSilfvergrenToFedMeal.reactions];
simSilfvergren.states         = [SimSilfvergrenToFedMeal.states];

% Declare fed state at first meal
simGly =   SimSilfvergrenFedMeal.variablevalues(1,ismember(SimSilfvergrenFedMeal.variables,'Glycogen_liver'));
if simGly < 250
    Cost_Model = Cost_Model + 5000;
end

simG   = simSilfvergren.variablevalues(:,ismember(simSilfvergren.variables,'G'))/18;
simG   = simG(round(Silfvergren2021_data.time(1:209),0));

% Create glucose weight used in parameter estimation
GlucoseWeight = Silfvergren2021_data.p1_glucoseCalibrated(1:209)*0.15;  % No SEM so a qualitative weight is added
GlucoseWeight(2:28)    = GlucoseWeight(2:28)*0.3;                       % Higher weight at regions of interest - fed state.
GlucoseWeight(2)       = GlucoseWeight(2)*0.4;                          % increased weight on initial value
GlucoseWeight(170:209) = GlucoseWeight(170:209)*0.3;                    % Higher weight at regions of interest - unfed state

CostG   = ((Silfvergren2021_data.p2_glucoseCalibrated(1:209) - simG(1:209)).^2)  ./(GlucoseWeight(1:209).^2); %No weight to initial datapoint

Cost_Model =  Cost_Model + sum(CostG,'omitnan');

end

