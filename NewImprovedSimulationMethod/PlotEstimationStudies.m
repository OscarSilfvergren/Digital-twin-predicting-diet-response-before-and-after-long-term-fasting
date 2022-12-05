%% pre-calculations
load('EstimationData_params');

AmountParametersOptimized = 58;
params = EstimationData_params;

CapSize            = 15;
LineWidthValue     = 15;
MealLineWidthValue = 15;
%% Declare individual params

Cost_Model = EstimationData_costfunction(Krssak2004Healthy_data,Rothman1991_data,Silfvergren2021_data,log(params),model,initialvalues, pNames,sNames,AmountParametersOptimized)

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
    SimKrssakAfterMeal             = model(375:1:1700,SimKrssak_meal.statevalues(end,:),params(1:65));             % Sim  after meal
    
 catch error
     Cost_Model = 1e60;
     return
 end

% Merge
simKrssak2004.time           = (1:1:1700)-360;
simKrssak2004.variablevalues = [SimKrssakToMeal.variablevalues(1:end-1,:); SimKrssak_meal.variablevalues(1:end-1,:); SimKrssakAfterMeal.variablevalues(1:end,:)];
simKrssak2004.reactionvalues = [SimKrssakToMeal.reactionvalues(1:end-1,:); SimKrssak_meal.reactionvalues(1:end-1,:); SimKrssakAfterMeal.reactionvalues(1:end,:)];
simKrssak2004.statevalues    = [SimKrssakToMeal.statevalues(1:end-1,:); SimKrssak_meal.statevalues(1:end-1,:); SimKrssakAfterMeal.statevalues(1:end,:)];
simKrssak2004.variables      = [SimKrssakToMeal.variables];
simKrssak2004.reactions      = [SimKrssakToMeal.reactions];
simKrssak2004.states         = [SimKrssakToMeal.states];

% --- Plot --- %


% Glycogen
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
grid on
set(gca,'xtick',[0,10,20],'ytick',[100,200,300],'FontSize', 70,'fontname','Arial')
plot(simKrssak2004.time/60, simKrssak2004.variablevalues(:,ismember(simKrssak2004.variables,'Glycogen_liver')),'b','LineWidth',LineWidthValue);
ylabel({'Hepatic glycogen' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
errorbar(Krssak2004Healthy_data.timeGlycogen(1:22)/60,Krssak2004Healthy_data.glycogen(1:22),Krssak2004Healthy_data.glycogenSEM(1:22),' k .','MarkerSize',1,'LineWidth',10,'CapSize',CapSize');
line([0 0.2], [150 150],'Color','k','LineWidth',MealLineWidthValue);
xlim([-0.5 20])
ylim([150 360])
hold off

% EGP
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
grid on
set(gca,'xtick',[0,5,10,15,20],'ytick',[0,0.5,1,1.5,2,2.5],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
plot(simKrssak2004.time/60,simKrssak2004.reactionvalues(:,ismember(simKrssak2004.reactions,'EGP')),'b','LineWidth',LineWidthValue);
ylabel({'EGP' ; '(mg/kg/min)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
errorbar(Krssak2004Healthy_data.timeEGP(1:20)/60,Krssak2004Healthy_data.EGP(1:20),Krssak2004Healthy_data.EGPSEM(1:20),' k .','MarkerSize',1,'LineWidth',10,'CapSize',CapSize');
line([0 0.2], [0.2 0.2],'Color','k','LineWidth',MealLineWidthValue);
xlim([-0.5 10])
ylim([0.2 1.85])
hold off

% G
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
grid on
a = gca;
set(gca,'xtick',[0,5,10,15,20],'ytick',[0,5,10,15,20],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
plot(simKrssak2004.time/60, simKrssak2004.variablevalues(:,ismember(simKrssak2004.variables,'G'))/18,'b','LineWidth',LineWidthValue)
ylabel({'Plasma glucose' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
errorbar(Krssak2004Healthy_data.timeGlucose(1:25)/60,Krssak2004Healthy_data.glucose(1:25)/18,Krssak2004Healthy_data.glucoseSEM(1:25)/18,' k .','MarkerSize',1,'LineWidth',10,'CapSize',CapSize');
line([0 0.2], [2 2],'Color','k','LineWidth',MealLineWidthValue);
xlim([-0.5 10])
ylim([2 11])
hold off

% I
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
grid on
a = gca;
set(gca,'xtick',[0,5,10,15,20],'ytick',[0,200,400,600,800,1000,1200],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
plot(simKrssak2004.time/60, simKrssak2004.variablevalues(:,ismember(simKrssak2004.variables,'I')),'b','LineWidth',LineWidthValue)
ylabel({'Plasma insulin' ; '(pM)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
errorbar(Krssak2004Healthy_data.timeInsulin(1:10)/60,Krssak2004Healthy_data.insulin(1:10),Krssak2004Healthy_data.insulinSEM(1:10),' k .','MarkerSize',1,'LineWidth',10,'CapSize',CapSize');
line([0 0.2], [0 0],'Color','k','LineWidth',MealLineWidthValue);
xlim([-0.5 10])
ylim([0 850])
hold off

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

% No need to merge as simulation is done in one step (no meals)
 
% --- Plot --- %

%Glycogen
figure('Name', "Glycogen", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
grid on
set(gca,'xtick',[0,30,60],'ytick',[0,250,500],'FontSize', 70,'FontSmoothing','on','fontname','Arial')
plot(SimRothman.time/60, SimRothman.variablevalues(:,ismember(SimRothman.variables,'Glycogen_liver')),'b','LineWidth',LineWidthValue);
errorbar(Rothman1991_data.timepALL(1:10)/60,Rothman1991_data.glycogenpALL(1:10),Rothman1991_data.glycogenpSEM(1:10),' k .','MarkerSize',1,'LineWidth',10,'CapSize',40');
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
ylabel({'Hepatic glycogen' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
xlim([0 62])
ylim([0 500])
hold off


%Gluconeogenesis
figure('Name', "Glycogen", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
grid on
set(gca,'xtick',[20,40,60],'ytick',[0,20,40,60,80,100],'FontSize', 70,'FontSmoothing','on','fontname','Arial')
Gluconeogenesis1    = sum(SimRothman.reactionvalues(1:1320,ismember(SimRothman.reactions,'Gluconeogenesis')));
Gluconeogenesis2    = sum(SimRothman.reactionvalues(1320:2160,ismember(SimRothman.reactions,'Gluconeogenesis')));
Gluconeogenesis3    = sum(SimRothman.reactionvalues(2160:3240,ismember(SimRothman.reactions,'Gluconeogenesis')));
EGP1                = sum(SimRothman.reactionvalues(1:1320,ismember(SimRothman.reactions,'EGP')));
EGP2                = sum(SimRothman.reactionvalues(1320:2160,ismember(SimRothman.reactions,'EGP')));
EGP3                = sum(SimRothman.reactionvalues(2160:3240,ismember(SimRothman.reactions,'EGP')));
Procentage1         = Gluconeogenesis1/EGP1;
Procentage2         = Gluconeogenesis2/EGP2;
Procentage3         = Gluconeogenesis3/EGP3;
Procentage          = [Procentage1, Procentage2, Procentage3];
temptime            = [1320, 2160, 3240];
plot(temptime/60, Procentage*100,'b','LineWidth',LineWidthValue)
errorbar(temptime/60,[0.64, 0.82, 0.96]*100,[0.05, 0.05, 0.01]*100,' k .','MarkerSize',1,'LineWidth',10,'CapSize',40');
xlim([20 60])
ylim([40 100])
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
ylabel({'Gluconeogenesis' ; 'contribution to EGP(%)'},'FontSmoothing','on','fontname','Arial');
hold off

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
paramsSilfvergrenMeal(ismember(pNames,'meal_solid'))     = 0;
paramsSilfvergrenMeal(ismember(pNames,'meal_liquid'))    = 1;

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
    SimSilfvergrenAfterUnfedMeal  = model(2884:1:3250,SimSilfvergrenUnfedMeal.statevalues(end,:),params(1:65));                   % sim after unfed meal
    
catch error
    Cost_Model = 1e60;
    return
end

% Merge
simSilfvergren.time           = 1:1:3250;
simSilfvergren.variablevalues = [SimSilfvergrenToFedMeal.variablevalues(1:end-1,:); SimSilfvergrenFedMeal.variablevalues(1:end-1,:); SimSilfvergrenToUnfedMeal.variablevalues(1:end-1,:); SimSilfvergrenUnfedMeal.variablevalues(1:end-1,:); SimSilfvergrenAfterUnfedMeal.variablevalues(1:end,:)];
simSilfvergren.reactionvalues = [SimSilfvergrenToFedMeal.reactionvalues(1:end-1,:); SimSilfvergrenFedMeal.reactionvalues(1:end-1,:); SimSilfvergrenToUnfedMeal.reactionvalues(1:end-1,:); SimSilfvergrenUnfedMeal.reactionvalues(1:end-1,:); SimSilfvergrenAfterUnfedMeal.reactionvalues(1:end,:)];
simSilfvergren.statevalues    = [SimSilfvergrenToFedMeal.statevalues              ; SimSilfvergrenFedMeal.statevalues(1:end-1,:)   ; SimSilfvergrenToUnfedMeal.statevalues(1:end-1,:)   ; SimSilfvergrenUnfedMeal.statevalues(1:end-1,:)   ; SimSilfvergrenAfterUnfedMeal.statevalues(1:end,:)];
simSilfvergren.variables      = [SimSilfvergrenToFedMeal.variables];
simSilfvergren.reactions      = [SimSilfvergrenToFedMeal.reactions];
simSilfvergren.states         = [SimSilfvergrenToFedMeal.states];


% --- Plot --- %

figure('Name', "Glucose fasting participant 1", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
grid on
a = gca;
set(a,'xtick',[0,1,2],'ytick',[0,2,3,4,5,6,7],'FontSize', 70,'fontname','Arial')
plot(Silfvergren2021_data.time(1:209)/1440,Silfvergren2021_data.p2_glucoseCalibrated(1:209),'k.','MarkerSize',70,'LineWidth',10)
plot(simSilfvergren.time/1440, simSilfvergren.variablevalues(:,ismember(simSilfvergren.variables,'G'))/18,'b','LineWidth',LineWidthValue);
line([29 49]/1440, [2 2],'Color','k','LineWidth',MealLineWidthValue);
line([2881 2901]/1440, [2 2],'Color','k','LineWidth',MealLineWidthValue);
ylabel({'Plasma glucose' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (days)",'FontSmoothing','on','fontname','Arial');
xlim([0 2.3])
ylim([2 6])
hold off
