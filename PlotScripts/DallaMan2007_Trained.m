% Old
modelName = 'DallamanOldData';
model = IQMmodel([modelName '.txt']);
IQMmakeMEXmodel(model);
model = str2func(modelName);
OptimizeDallaManOld_params(41:43) = [0,60,75000];
simOld = model(time,[],OptimizeDallaManOld_params);

% New
load('DallaMan2007_param');
modelName = 'ModelSimpleMetabolismSteadyState';
model = IQMmodel([modelName '.txt']);
IQMmakeMEXmodel(model);
model = str2func(modelName);

%Sort for cost
clear optimizedParamTemp2
[row column] = size(DallaMan2007_Params);
k = 1;
for i = 1:row
    optimizedParamTemp  = DallaMan2007_Params(i,1:(column-1));
    cost                = DallaMan2007_costfunctionHealthy(DallaMan2007_data,time,log(optimizedParamTemp),modelName, meal_information,body_information);
    
    if cost < chi2inv(0.95,190)
        optimizedParamTemp2(k,column)     = cost;
        optimizedParamTemp2(k,1:column-1) = optimizedParamTemp;
        k = k+1;
    end
end

DallaMan2007_Params = sortrows(optimizedParamTemp2,column);
[row column] = size(DallaMan2007_Params);

disp(' ')
if DallaMan2007_Params(1,column) < chi2inv(0.95,126)
    disp('----- Best fit to Dalla Man data is below treshold -----')
end

fprintf('BestFit DallaManData: %.2f, Statistical Limit: %.2f (dgf = %i)', DallaMan2007_Params(1,column), chi2inv(0.95,126), 126)
disp(' ')

%%

simNew = model(time,[],DallaMan2007_Params(1,1:column-1));

% Glucose
figure('Name', "Compare Glucose", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
set(gca,'xtick',[0,4,8],'ytick',[4,8,12],'FontSize', 70,'FontSmoothing','on','fontname','Arial')
plot(time/60-127.5,simOld.reactionvalues(:,ismember(simOld.reactions,'G'))/18,'Color',colorOrange,'LineWidth',LineWidthValue)
plot(time/60-127.5,simNew.variablevalues(:,ismember(simNew.variables,'G'))/18,'Color',colorGrey,'LineWidth',LineWidthValue)
xlabel("Time (h)",'FontSmoothing','on');
ylabel({'Plasma glucose' ; '(mM)'},'FontSmoothing','on');
errorbar((DallaMan2007_data.Glucose_time(1:21)+7680)/60-127.5,DallaMan2007_data.Glucose(1:21)/18,DallaMan2007_data.GlucoseSEM(1:21)/18,' k .','MarkerSize',8,'LineWidth',4);
line([0.6 0.8], [4 4],'Color','k','LineWidth',8);
xlim([0 8])
ylim([4 12])
hold off


% Insulin
figure('Name', "Compare Insulin", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
set(gca,'xtick',[0,4,8],'ytick',[0,250,500],'FontSize', 70,'FontSmoothing','on','fontname','Arial')
plot(time/60-127.5,simOld.reactionvalues(:,ismember(simOld.reactions,'I')),'Color',colorOrange,'LineWidth',LineWidthValue)
plot(time/60-127.5,simNew.variablevalues(:,ismember(simNew.variables,'I')),'Color',colorGrey,'LineWidth',LineWidthValue)
xlabel("Time (h)",'FontSmoothing','on');
ylabel({'Plasma insulin' ; '(pM)'},'FontSmoothing','on');
line([0.6 0.8], [0 0],'Color','k','LineWidth',8);
errorbar((DallaMan2007_data.Insulin_time(1:21)+7680)/60-127.5,DallaMan2007_data.Insulin(1:21),DallaMan2007_data.InsulinSEM(1:21), " k . ",'MarkerSize',8,'LineWidth',4);
xlim([0 8])
ylim([0 600])
hold off

% EGP
figure('Name', "Compare EGP", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
set(gca,'xtick',[0,4,8],'ytick',[0,1.5,3],'FontSize', 70,'FontSmoothing','on','fontname','Arial')
plot(time/60-127.5,simOld.reactionvalues(:,ismember(simOld.reactions,'EGP')),'Color',colorOrange,'LineWidth',LineWidthValue)
plot(time/60-127.5,simNew.reactionvalues(:,ismember(simNew.reactions,'EGP')),'Color',colorGrey,'LineWidth',LineWidthValue)
xlabel("Time (h)",'FontSmoothing','on');
ylabel({'EGP' ; '(mg/kg/min)'},'FontSmoothing','on');
line([0.6 0.8], [0 0],'Color','k','LineWidth',8);
errorbar((DallaMan2007_data.EGP_time(1:22)+7680)/60-127.5,DallaMan2007_data.EGP(1:22),DallaMan2007_data.EGPSEM(1:22), " k . ",'MarkerSize',8,'LineWidth',4);
xlim([0 8])
ylim([0 3])
hold off

% Ra
figure('Name', "Compare Ra", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
set(gca,'xtick',[0,4,8],'ytick',[0,7,14],'FontSize', 70,'FontSmoothing','on','fontname','Arial')
plot(time/60-127.5,simOld.reactionvalues(:,ismember(simOld.reactions,'Ra')),'Color',colorOrange,'LineWidth',LineWidthValue)
plot(time/60-127.5,simNew.reactionvalues(:,ismember(simNew.reactions,'Ra')),'Color',colorGrey,'LineWidth',LineWidthValue)
xlabel("Time (h)",'FontSmoothing','on');
ylabel({'Rate of appearance' ; '(mg/kg/min)'},'FontSmoothing','on');
line([0.6 0.8], [0 0],'Color','k','LineWidth',8);
errorbar((DallaMan2007_data.Ra_time(1:20)+7680)/60-127.5,DallaMan2007_data.Ra(1:20),DallaMan2007_data.RaSEM(1:20), " k . ",'MarkerSize',8,'LineWidth',4);
xlim([0 8])
ylim([0 14])
hold off

% U_idii
figure('Name', "Compare glucose uptake", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
set(gca,'xtick',[0,4,8],'ytick',[0,4,8],'FontSize', 70,'FontSmoothing','on','fontname','Arial')
plot(time/60-127.5,simOld.reactionvalues(:,ismember(simOld.reactions,'U')),'Color',colorOrange,'LineWidth',LineWidthValue)
plot(time/60-127.5,simNew.reactionvalues(:,ismember(simNew.reactions,'U_idii')),'Color',colorGrey,'LineWidth',LineWidthValue)
xlabel("Time (h)",'FontSmoothing','on');
ylabel({'Glucose utilization' ; '(mg/kg/min)'},'FontSmoothing','on');
line([0.6 0.8], [0 0],'Color','k','LineWidth',8);
errorbar((DallaMan2007_data.U_time(1:21)+7680)/60-127.5,DallaMan2007_data.U(1:21),DallaMan2007_data.USEM(1:21), " k . ",'MarkerSize',8,'LineWidth',4);
xlim([0 8])
hold off

% S
figure('Name', "Compare Insulin S", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
set(gca,'xtick',[0,4,8],'ytick',[0,7,14],'FontSize', 70,'FontSmoothing','on','fontname','Arial')
plot(time/60-127.5,simOld.reactionvalues(:,ismember(simOld.reactions,'S')),'Color',colorOrange,'LineWidth',LineWidthValue)
plot(time/60-127.5,simNew.reactionvalues(:,ismember(simNew.reactions,'S')),'Color',colorGrey,'LineWidth',LineWidthValue)
xlabel("Time (h)",'FontSmoothing','on');
ylabel({'Insulin secration' ; '(pmol/kg/min)'},'FontSmoothing','on');
line([0.6 0.8], [0 0],'Color','k','LineWidth',8);
errorbar(((DallaMan2007_data.S_time(1:21)+7680))/60-127.5,DallaMan2007_data.S(1:21),DallaMan2007_data.SSEM(1:21), " k . ",'MarkerSize',8,'LineWidth',4);
xlim([0 8])
ylim([0 14])
hold off

%% Model Uncertainty unmeasured variables
time2 = [time, fliplr(time)];

[row column] = size(DallaMan2007_Params);
DallaMan2007_Params = sortrows(DallaMan2007_Params,column);


for i = 1:row
    optimizedParamTemp = DallaMan2007_Params(i,1:(column-1));
    sim = model(time,[],optimizedParamTemp);
    
    if i == 1
        simBest = sim;
    end
    % S
    if i == 1
        maxGlycogen_liver1 = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
        minGlycogen_liver1 = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
        maxGlycogen_liver2 = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
        minGlycogen_liver2 = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
    else
        maxGlycogen_liver1 = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
        minGlycogen_liver1 = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
        maxGlycogen_liver2 = max(maxGlycogen_liver2,maxGlycogen_liver1);
        minGlycogen_liver2 = min(minGlycogen_liver2,minGlycogen_liver1);
    end
    
    % S
    if i == 1
        maxGluconeogenesis1 = sim.reactionvalues(:,ismember(sim.reactions,'Gluconeogenesis'));
        minGluconeogenesis1 = sim.reactionvalues(:,ismember(sim.reactions,'Gluconeogenesis'));
        maxGluconeogenesis2 = sim.reactionvalues(:,ismember(sim.reactions,'Gluconeogenesis'));
        minGluconeogenesis2 = sim.reactionvalues(:,ismember(sim.reactions,'Gluconeogenesis'));
    else
        maxGluconeogenesis1 = sim.reactionvalues(:,ismember(sim.reactions,'Gluconeogenesis'));
        minGluconeogenesis1 = sim.reactionvalues(:,ismember(sim.reactions,'Gluconeogenesis'));
        maxGluconeogenesis2 = max(maxGluconeogenesis2,maxGluconeogenesis1);
        minGluconeogenesis2 = min(minGluconeogenesis2,minGluconeogenesis1);
    end
    
    % S
    if i == 1
        maxEGP_Kidneys1 = sim.reactionvalues(:,ismember(sim.reactions,'EGP_Kidneys'));
        minEGP_Kidneys1 = sim.reactionvalues(:,ismember(sim.reactions,'EGP_Kidneys'));
        maxEGP_Kidneys2 = sim.reactionvalues(:,ismember(sim.reactions,'EGP_Kidneys'));
        minEGP_Kidneys2 = sim.reactionvalues(:,ismember(sim.reactions,'EGP_Kidneys'));
    else
        maxEGP_Kidneys1 = sim.reactionvalues(:,ismember(sim.reactions,'EGP_Kidneys'));
        minEGP_Kidneys1 = sim.reactionvalues(:,ismember(sim.reactions,'EGP_Kidneys'));
        maxEGP_Kidneys2 = max(maxEGP_Kidneys2,maxEGP_Kidneys1);
        minEGP_Kidneys2 = min(minEGP_Kidneys2,minEGP_Kidneys1);
    end
    
end
%%
% Glycogen
figure('Name', "Predicted Gly", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
set(gca,'xtick',[0,6,12],'ytick',[100,200,300,400],'FontSize', 70,'FontSmoothing','on','fontname','Arial')
y   = [maxGlycogen_liver2', fliplr(minGlycogen_liver2')];
line([0.6 0.8], [100 100],'Color','k','LineWidth',8);
fill(time2/60-127.5,y,colorGrey,'FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/60-127.5,simBest.variablevalues(:,ismember(simBest.variables,'Glycogen_liver')),'Color',colorGrey,'LineWidth',LineWidthValue)
xlabel("Time (h)",'FontSmoothing','on');
ylabel({'Hepatic glycogen' ; '(mM)'},'FontSmoothing','on');
xlim([0 12])
ylim([100 400])
hold off

% Gluconeogenesis
figure('Name', "Predicted GNG", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
set(gca,'xtick',[0,6,12],'ytick',[0,0.2,0.4,0.6,0.8],'FontSize', 70,'FontSmoothing','on','fontname','Arial')
y   = [maxGluconeogenesis2', fliplr(minGluconeogenesis2')];
fill(time2/60-127.5,y,colorGrey,'FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/60-127.5,simBest.reactionvalues(:,ismember(simBest.reactions,'Gluconeogenesis')),'Color',colorGrey,'LineWidth',LineWidthValue)
line([0.6 0.8], [0 0],'Color','k','LineWidth',8);
xlabel("Time (h)",'FontSmoothing','on');
ylabel({'Gluconeogenesis' ; '(mg/kg/min)'},'FontSmoothing','on');
xlim([0 12])
ylim([0 0.8])
hold off

% Kidneys
figure('Name', "Predicted renal EGP", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
set(gca,'xtick',[0,6,12],'ytick',[0.2,0.25,0.3],'FontSize', 70,'FontSmoothing','on','fontname','Arial')
y   = [maxEGP_Kidneys2', fliplr(minEGP_Kidneys2')];
fill(time2/60-127.5,y,colorGrey,'FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/60-127.5,simBest.reactionvalues(:,ismember(simBest.reactions,'EGP_Kidneys')),'Color',colorGrey,'LineWidth',LineWidthValue)
line([0.6 0.8], [0.18 0.18],'Color','k','LineWidth',8);
xlabel("Time (h)",'FontSmoothing','on');
ylabel({'Renal EGP' ; '(mg/kg/min)'},'FontSmoothing','on');
xlim([0 12])
ylim([0.18 0.32])
hold off

%% Model uncertainty ratios
time2 = [time, fliplr(time)];

clear ratioK
clear ratioL
clear ratio_f
clear ratio_m
clear ratio_b
clear ratio_l
clear ratioLiver
clear ratioBlood

[row column] = size(DallaMan2007_Params);
DallaMan2007_Params = sortrows(DallaMan2007_Params,column);

for i = 1:row
    optimizedParamTemp = DallaMan2007_Params(i,1:(column-1));
    sim = model(time,[],optimizedParamTemp);
    
    EGP_Kidneys = sum(sim.reactionvalues(7600:7960,ismember(sim.reactions,'EGP_Kidneys')));
    EGP_Liver   = sum(sim.reactionvalues(7600:7960,ismember(sim.reactions,'EGP_Liver')));
    EGP         = sum(sim.reactionvalues(7600:7960,ismember(sim.reactions,'EGP')));
    
    U_idl = sum(sim.reactionvalues(7240:7600,ismember(sim.reactions,'U_il')));
    U_idm = sum(sim.reactionvalues(7240:7600,ismember(sim.reactions,'U_idm')));
    U_idf = sum(sim.reactionvalues(7240:7600,ismember(sim.reactions,'U_idf')));
    U_ii  = sum(sim.reactionvalues(7240:7600,ismember(sim.reactions,'U_ii_print')));
    U     = sum(sim.reactionvalues(7240:7600,ismember(sim.reactions,'U')));
    
    InsulinDegradation_Blood     = sum(sim.reactionvalues(7680:7960,ismember(sim.reactions,'InsulinDegradation_Blood')));
    InsulinDegradation_Liver     = sum(sim.reactionvalues(7680:7960,ismember(sim.reactions,'InsulinDegradation_Liver')));
    
    ratioK(i)       = EGP_Kidneys/EGP;
    ratioL(i)       = EGP_Liver/EGP;
    
    ratio_f(i) = U_idf/U;
    ratio_m(i) = U_idm/U;
    ratio_b(i) = U_ii/U;
    ratio_l(i) = U_idl/U;
    
    ratioLiver(i)  = InsulinDegradation_Liver/(InsulinDegradation_Blood + InsulinDegradation_Liver); % 0.5 till 0.8 "percentage varies widely"
    ratioBlood(i)  = InsulinDegradation_Blood/(InsulinDegradation_Blood + InsulinDegradation_Liver); % 0.5 till 0.8 "percentage varies widely"
    
end
%%
% EGP vs kidneys
figure('Name', "EGP vs kidney", 'units', 'normalized', 'outerposition', [0 0 1 1])

% X = categorical({'Sim Kidneys','Data Kidneys','Sim Liver','Data Liver'});
% X = reordercats(X,{'Sim Kidneys','Data Kidneys','Sim Liver','Data Liver'});

X = categorical({'Sim 1','Data 1','Sim 2','Data 2'});
X = reordercats(X,{'Sim 1','Data 1','Sim 2','Data 2'});


KidneysBestMinMax = BestMinMax(ratioK,0.2);
LiverBestMinMax   = BestMinMax(ratioL,0.8);

hold on
set(gca,'ytick',[0,50,100],'FontSize', 70,'FontSmoothing','on','fontname','Arial')
bar(X,[KidneysBestMinMax(1);20;LiverBestMinMax(1);80])
errorbar(X,[KidneysBestMinMax(1);20;LiverBestMinMax(1);80],[KidneysBestMinMax(2),10,LiverBestMinMax(2),10], " k . ",'MarkerSize',8,'LineWidth',4);
ylabel({'Contribution to EGP' ; '(%)'},'FontSmoothing','on');
ylim([0 100])
hold off
%%
% Fasting; U
figure('Name', "glucose uptake organ", 'units', 'normalized', 'outerposition', [0 0 1 1])

 X = categorical({'S Fat','D Fat','S Muscle','D Muscle','S Brain','D Brain','S Liver','D Liver'});
 X = reordercats(X,{'S Fat','D Fat','S Muscle','D Muscle','S Brain','D Brain','S Liver','D Liver'});

FatBestMinMax     = BestMinMax(ratio_f,17.5/100);
MuscleBestMinMax  = BestMinMax(ratio_m,22.5/100);
BrainBestMinMax   = BestMinMax(ratio_b,40/100);
LiverBestMinMax   = BestMinMax(ratio_l,20/100);

hold on
set(gca,'ytick',[0,50,100],'FontSize', 70,'FontSmoothing','on','fontname','Arial')
bar(X,[[FatBestMinMax(1);17.5;MuscleBestMinMax(1);22.5;BrainBestMinMax(1);40;LiverBestMinMax(1);20]])
errorbar(X,[FatBestMinMax(1);17.5;MuscleBestMinMax(1);22.5;BrainBestMinMax(1);40;LiverBestMinMax(1);20],[[FatBestMinMax(2);7.5;MuscleBestMinMax(2);7.5;BrainBestMinMax(2);7.5;LiverBestMinMax(2);7.5]], " k . ",'MarkerSize',8,'LineWidth',4);
ylabel({'Organ specific glucose' ; 'uptake during fasting (%)'},'FontSmoothing','on');
ylim([0 100])
hold off

%% Insulin clearance
figure('Name', "Insulin clearance uptake organ", 'units', 'normalized', 'outerposition', [0 0 1 1])

X = categorical({'S Liver','D Liver','S Non-Liver','D Non-Liver'});
X = reordercats(X,{'S Liver','D Liver','S Non-Liver','D Non-Liver'});

LiverBestMinMax     = BestMinMax(ratioLiver,65/100);
BloodBestMinMax  = BestMinMax(ratioBlood,25/100);

hold on
set(gca,'ytick',[0,50,100],'FontSize', 70,'FontSmoothing','on','fontname','Arial')
bar(X,[LiverBestMinMax(1); 65 ; BloodBestMinMax(1);25]);
errorbar(X,[LiverBestMinMax(1); 65 ; BloodBestMinMax(1);25],[LiverBestMinMax(2); 10 ; BloodBestMinMax(2);10], " k . ",'MarkerSize',8,'LineWidth',4);
ylabel({'Insulin clearance' ; '(%)'},'FontSmoothing','on');
ylim([0 100])
hold off