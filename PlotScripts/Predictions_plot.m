%% pre-calc
load('ModelStartGuess');
time = linspace(0,28480,28481);
[row column] = size(ModelStartGuess);
params = ModelStartGuess;
resett = ModelStartGuess;
startTime = 12960;
endTime   = 22920;
diff = endTime - startTime;
disp(' ')
disp('--- Prediction of diets and different situations ---')
disp(' ')

%% Intermediate fasting, 2 meals; 2000kcal:45carb/27.5protein/27.5fett
modelName = 'ModelSimpleMetabolismIntermediatefasting';
model = IQMmodel([modelName '.txt']);
IQMmakeMEXmodel(model);
model = str2func(modelName);

params = resett;
params(110:113) = [30, 0, 112.5, 68.75];
simfast = model(time,[],params);

G               = simfast.variablevalues((startTime:endTime),ismember(simfast.variables,'G'))/18;
Gly             = simfast.variablevalues((startTime:endTime),ismember(simfast.variables,'Glycogen_liver'));
EGP             = simfast.reactionvalues((startTime:endTime),ismember(simfast.reactions,'EGP'));
I               = simfast.variablevalues((startTime:endTime),ismember(simfast.variables,'I'));
Gluconeogenesis = simfast.reactionvalues((startTime:endTime),ismember(simfast.reactions,'Gluconeogenesis'));

Gsave(1)               = mean(G);
EGPsave(1)             = mean(EGP);
Isave(1)               = mean(I);
Gluconeogenesissave(1) = mean(Gluconeogenesis);
Glysave(1)             = mean(Gly);

%% 5:2 Vanlig dag: 2560 kcal:45carb/27.5protein/27.5fett  lowdag: 600kcal:45carb/27.5protein/27.5fett
modelName = 'ModelSimpleMetabolism5_2';
model = IQMmodel([modelName '.txt']);
IQMmakeMEXmodel(model);
model = str2func(modelName);

params = resett;
params(77:80)  = [15, 0, 33.75, 20.625];  %Low
params(110:113) = [30, 0, 96, 58.6];   %High
sim52 = model(time,[],params);

G               = sim52.variablevalues((startTime:endTime),ismember(sim52.variables,'G'))/18;
Gly             = sim52.variablevalues((startTime:endTime),ismember(sim52.variables,'Glycogen_liver'));
EGP             = sim52.reactionvalues((startTime:endTime),ismember(sim52.reactions,'EGP'));
I               = sim52.variablevalues((startTime:endTime),ismember(sim52.variables,'I'));
Gluconeogenesis = sim52.reactionvalues((startTime:endTime),ismember(sim52.reactions,'Gluconeogenesis'));

Gsave(2)               = mean(G);
EGPsave(2)             = mean(EGP);
Isave(2)               = mean(I);
Gluconeogenesissave(2) = mean(Gluconeogenesis);
Glysave(2)             = mean(Gly);

%% SFM, 4 meals; 2000kcal:45carb/27.5protein/27.5fett
modelName = 'ModelSimpleMetabolismSFM';
model = IQMmodel([modelName '.txt']);
IQMmakeMEXmodel(model);
model = str2func(modelName);

params = resett;
params(110:113) = [15, 0, 56.25, 34.375];
simSFM = model(time,[],params);

G               = simSFM.variablevalues((startTime:endTime),ismember(simSFM.variables,'G'))/18;
Gly             = simSFM.variablevalues((startTime:endTime),ismember(simSFM.variables,'Glycogen_liver'));
EGP             = simSFM.reactionvalues((startTime:endTime),ismember(simSFM.reactions,'EGP'));
I               = simSFM.variablevalues((startTime:endTime),ismember(simSFM.variables,'I'));
Gluconeogenesis = simSFM.reactionvalues((startTime:endTime),ismember(simSFM.reactions,'Gluconeogenesis'));

Gsave(3)               = mean(G);
EGPsave(3)             = mean(EGP);
Isave(3)               = mean(I);
Gluconeogenesissave(3) = mean(Gluconeogenesis);
Glysave(3)             = mean(Gly);

%% LCHF  2000kcal:30carb/35protein/35fett
modelName = 'ModelSimpleMetabolismPredictionDiet';
model = IQMmodel([modelName '.txt']);
IQMmakeMEXmodel(model);
model = str2func(modelName);

params = resett;
params(110:113) = [20, 0, 50, 58];
simLCHF = model(time,[],params);

G               = simLCHF.variablevalues((startTime:endTime),ismember(simLCHF.variables,'G'))/18;
Gly             = simLCHF.variablevalues((startTime:endTime),ismember(simLCHF.variables,'Glycogen_liver'));
EGP             = simLCHF.reactionvalues((startTime:endTime),ismember(simLCHF.reactions,'EGP'));
I               = simLCHF.variablevalues((startTime:endTime),ismember(simLCHF.variables,'I'));
Gluconeogenesis = simLCHF.reactionvalues((startTime:endTime),ismember(simLCHF.reactions,'Gluconeogenesis'));

Gsave(4)               = mean(G);
EGPsave(4)             = mean(EGP);
Isave(4)               = mean(I);
Gluconeogenesissave(4) = mean(Gluconeogenesis);
Glysave(4)             = mean(Gly);

%% HCLF 2000kcal:60carb/25protein/15fett
modelName = 'ModelSimpleMetabolismPredictionDiet';
model = IQMmodel([modelName '.txt']);
IQMmakeMEXmodel(model);
model = str2func(modelName);

params = resett;
params(110:113) = [60, 0, 100, 41];
simHCLF = model(time,[],params);

G               = simHCLF.variablevalues((startTime:endTime),ismember(simHCLF.variables,'G'))/18;
Gly             = simHCLF.variablevalues((startTime:endTime),ismember(simHCLF.variables,'Glycogen_liver'));
EGP             = simHCLF.reactionvalues((startTime:endTime),ismember(simHCLF.reactions,'EGP'));
I               = simHCLF.variablevalues((startTime:endTime),ismember(simHCLF.variables,'I'));
Gluconeogenesis = simHCLF.reactionvalues((startTime:endTime),ismember(simHCLF.reactions,'Gluconeogenesis'));

Gsave(5)               = mean(G);
EGPsave(5)             = mean(EGP);
Isave(5)               = mean(I);
Gluconeogenesissave(5) = mean(Gluconeogenesis);
Glysave(5)             = mean(Gly);

%% Eating slow  2000kcal:50carb/25protein/25fett
modelName = 'ModelSimpleMetabolismPredictionDiet';
model = IQMmodel([modelName '.txt']);
IQMmakeMEXmodel(model);
model = str2func(modelName);

params = resett;
params(110:113) = [40, 0, 83.3, 41.6];
simSlow = model(time,[],params);

G               = simSlow.variablevalues((startTime:endTime),ismember(simSlow.variables,'G'))/18;
Gly             = simSlow.variablevalues((startTime:endTime),ismember(simSlow.variables,'Glycogen_liver'));
EGP             = simSlow.reactionvalues((startTime:endTime),ismember(simSlow.reactions,'EGP'));
I               = simSlow.variablevalues((startTime:endTime),ismember(simSlow.variables,'I'));
Gluconeogenesis = simSlow.reactionvalues((startTime:endTime),ismember(simSlow.reactions,'Gluconeogenesis'));

Gsave(6)               = mean(G);
EGPsave(6)             = mean(EGP);
Isave(6)               = mean(I);
Gluconeogenesissave(6) = mean(Gluconeogenesis);
Glysave(6)             = mean(Gly);

%% Eating fast  2000kcal:30carb/35protein/35fett
modelName = 'ModelSimpleMetabolismPredictionDiet';
model = IQMmodel([modelName '.txt']);
IQMmakeMEXmodel(model);
model = str2func(modelName);

params = resett;
params(110:113) = [10, 0, 83.3, 41.6];
simFast = model(time,[],params);

G               = simFast.variablevalues((startTime:endTime),ismember(simFast.variables,'G'))/18;
Gly             = simFast.variablevalues((startTime:endTime),ismember(simFast.variables,'Glycogen_liver'));
EGP             = simFast.reactionvalues((startTime:endTime),ismember(simFast.reactions,'EGP'));
I               = simFast.variablevalues((startTime:endTime),ismember(simFast.variables,'I'));
Gluconeogenesis = simFast.reactionvalues((startTime:endTime),ismember(simFast.reactions,'Gluconeogenesis'));

Gsave(7)               = mean(G);
EGPsave(7)             = mean(EGP);
Isave(7)               = mean(I);
Gluconeogenesissave(7) = mean(Gluconeogenesis);
Glysave(7)             = mean(Gly);

%% Big person  2000kcal:50carb/25protein/25fett
modelName = 'ModelSimpleMetabolismPredictionDiet';
model = IQMmodel([modelName '.txt']);
IQMmakeMEXmodel(model);
model = str2func(modelName);

params = resett;
params(102:105) = [0 , 1, 180, 100 ];  % female, male, height, weight
params(110:113) = [30, 0, 83.3, 41.6];
simBig = model(time,[],params);

G               = simBig.variablevalues((startTime:endTime),ismember(simBig.variables,'G'))/18;
Gly             = simBig.variablevalues((startTime:endTime),ismember(simBig.variables,'Glycogen_liver'));
EGP             = simBig.reactionvalues((startTime:endTime),ismember(simBig.reactions,'EGP'));
I               = simBig.variablevalues((startTime:endTime),ismember(simBig.variables,'I'));
Gluconeogenesis = simBig.reactionvalues((startTime:endTime),ismember(simBig.reactions,'Gluconeogenesis'));

Gsave(8)               = mean(G);
EGPsave(8)             = mean(EGP);
Isave(8)               = mean(I);
Gluconeogenesissave(8) = mean(Gluconeogenesis);
Glysave(8)             = mean(Gly);

%% Small person  2000kcal:30carb/35protein/35fett
modelName = 'ModelSimpleMetabolismPredictionDiet';
model = IQMmodel([modelName '.txt']);
IQMmakeMEXmodel(model);
model = str2func(modelName);

params = resett;
params(102:105) = [0 , 1, 180, 70 ];  % female, male, height, weight
params(110:113) = [30, 0, 83.3, 41.6];
simSmall = model(time,[],params);

G               = simSmall.variablevalues((startTime:endTime),ismember(simSmall.variables,'G'))/18;
Gly             = simSmall.variablevalues((startTime:endTime),ismember(simSmall.variables,'Glycogen_liver'));
EGP             = simSmall.reactionvalues((startTime:endTime),ismember(simSmall.reactions,'EGP'));
I               = simSmall.variablevalues((startTime:endTime),ismember(simSmall.variables,'I'));
Gluconeogenesis = simSmall.reactionvalues((startTime:endTime),ismember(simSmall.reactions,'Gluconeogenesis'));

Gsave(9)               = mean(G);
EGPsave(9)             = mean(EGP);
Isave(9)               = mean(I);
Gluconeogenesissave(9) = mean(Gluconeogenesis);
Glysave(9)             = mean(Gly);


%% Plot diets Gly
%Gly
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,1,2,3],'ytick',[0,200,400],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
plot((time-12960)/1440,simfast.variablevalues(:,ismember(simfast.variables,'Glycogen_liver')),'LineWidth',LineWidthValue)
ylabel({'Hepatic glycogen' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (days)",'FontSmoothing','on','fontname','Arial');
title('Intermittent fasting','FontSmoothing','on','fontname','Arial');
x1=[0.44,0.48];
x2=[0.60, 0.64];
line(x1, [0 0],'Color','k','LineWidth',12);
line(x2, [0 0],'Color','k','LineWidth',12);
line(x1+1, [0 0],'Color','k','LineWidth',12);
line(x2+1, [0 0],'Color','k','LineWidth',12);
line(x1+2, [0 0],'Color','k','LineWidth',12);
line(x2+2, [0 0],'Color','k','LineWidth',12);
xlim([0, 3])
ylim([0,400])
hold off


%Gly
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,1,2,3],'ytick',[0,200,400],'FontSize', 70,'fontname','Arial')
plot((time-12960)/1440,sim52.variablevalues(:,ismember(sim52.variables,'Glycogen_liver')),'LineWidth',LineWidthValue)
ylabel({'Hepatic glycogen' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (days)",'FontSmoothing','on','fontname','Arial');
title('5:2','FontSmoothing','on','fontname','Arial');
x1=[0.23,0.27];
x2=[0.47, 0.51];
x3=[0.78, 0.82];
line(x1, [0 0],'Color','k','LineWidth',12);
line(x2, [0 0],'Color','k','LineWidth',12);
line(x3, [0 0],'Color','k','LineWidth',12);
line(x1+1, [0 0],'Color','k','LineWidth',12);
line(x2+1, [0 0],'Color','k','LineWidth',12);
line(x3+1, [0 0],'Color','k','LineWidth',12);
line(x1+2, [0 0],'Color','k','LineWidth',12);
line(x3+2, [0 0],'Color','k','LineWidth',12);
xlim([0, 3])
ylim([0,400])
hold off


%Gly
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,1,2,3],'ytick',[0,200,400],'FontSize', 70,'fontname','Arial')
plot((time-12960)/1440,simSFM.variablevalues(:,ismember(simSFM.variables,'Glycogen_liver')),'LineWidth',LineWidthValue)
ylabel({'Hepatic glycogen' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (days)",'FontSmoothing','on','fontname','Arial');
title('Small frequent meals','FontSmoothing','on','fontname','Arial');
x1=[0.34,0.38];
x2=[0.50, 0.54];
x3=[0.67, 0.71];
x4=[0.83, 0.87];
line(x1, [0 0],'Color','k','LineWidth',12);
line(x2, [0 0],'Color','k','LineWidth',12);
line(x3, [0 0],'Color','k','LineWidth',12);
line(x4, [0 0],'Color','k','LineWidth',12);
line(x1+1, [0 0],'Color','k','LineWidth',12);
line(x2+1, [0 0],'Color','k','LineWidth',12);
line(x3+1, [0 0],'Color','k','LineWidth',12);
line(x4+1, [0 0],'Color','k','LineWidth',12);
line(x1+2, [0 0],'Color','k','LineWidth',12);
line(x2+2, [0 0],'Color','k','LineWidth',12);
line(x3+2, [0 0],'Color','k','LineWidth',12);
line(x4+2, [0 0],'Color','k','LineWidth',12);
xlim([0, 3])
ylim([0,400])
hold off


% Gly
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,1,2,3],'ytick',[0,200,400],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
plot((time-12960)/1440,simLCHF.variablevalues(:,ismember(simLCHF.variables,'Glycogen_liver')),'LineWidth',LineWidthValue)
ylabel({'Hepatic glycogen' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (days)",'FontSmoothing','on','fontname','Arial');
title('Low carbs high fat','FontSmoothing','on','fontname','Arial');
x1=[0.23,0.27];
x2=[0.47, 0.51];
x3=[0.78, 0.82];
line(x1, [0 0],'Color','k','LineWidth',12);
line(x2, [0 0],'Color','k','LineWidth',12);
line(x3, [0 0],'Color','k','LineWidth',12);
line(x1+1, [0 0],'Color','k','LineWidth',12);
line(x2+1, [0 0],'Color','k','LineWidth',12);
line(x3+1, [0 0],'Color','k','LineWidth',12);
line(x1+2, [0 0],'Color','k','LineWidth',12);
line(x2+2, [0 0],'Color','k','LineWidth',12);
line(x3+2, [0 0],'Color','k','LineWidth',12);
xlim([0, 3])
ylim([0,400])
hold off


% Gly
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,1,2,3],'ytick',[0,200,400],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
plot((time-12960)/1440,simHCLF.variablevalues(:,ismember(simHCLF.variables,'Glycogen_liver')),'LineWidth',LineWidthValue)
ylabel({'Hepatic glycogen' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (days)",'FontSmoothing','on','fontname','Arial');
title('High carbs low fat','FontSmoothing','on','fontname','Arial');
x1=[0.23,0.27];
x2=[0.47, 0.51];
x3=[0.78, 0.82];
line(x1, [0 0],'Color','k','LineWidth',12);
line(x2, [0 0],'Color','k','LineWidth',12);
line(x3, [0 0],'Color','k','LineWidth',12);
line(x1+1, [0 0],'Color','k','LineWidth',12);
line(x2+1, [0 0],'Color','k','LineWidth',12);
line(x3+1, [0 0],'Color','k','LineWidth',12);
line(x1+2, [0 0],'Color','k','LineWidth',12);
line(x2+2, [0 0],'Color','k','LineWidth',12);
line(x3+2, [0 0],'Color','k','LineWidth',12);
xlim([0, 3])
ylim([0,400])
hold off

%% Plot bars diets

X = categorical({'IF','5:2','SFM','LCHF','HCLF'});
X = reordercats(X,{'IF','5:2','SFM','LCHF','HCLF'});

% Glucose AUC
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'ytick',[0,3,6],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
bar(X,[Gsave(1);0;0;0;0])
bar(X,[0;Gsave(2);0;0;0])
bar(X,[0;0;Gsave(3);0;0])
bar(X,[0;0;0;Gsave(4);0])
bar(X,[0;0;0;0;Gsave(5)])
ylabel({'Mean plasma glucose' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
% ylim([3,5.7])
hold off
%%

% I
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'ytick',[0,150,300],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
bar(X,[Isave(1);0;0;0;0])
bar(X,[0;Isave(2);0;0;0])
bar(X,[0;0;Isave(3);0;0])
bar(X,[0;0;0;Isave(4);0])
bar(X,[0;0;0;0;Isave(5)])
ylabel({'Mean plasma insulin' ; '(pM)'},'FontSmoothing','on','fontname','Arial');
% ylim([80,240])
hold off
%%

% Gly
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'ytick',[0,150,300],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
bar(X,[Glysave(1);0;0;0;0])
bar(X,[0;Glysave(2);0;0;0])
bar(X,[0;0;Glysave(3);0;0])
bar(X,[0;0;0;Glysave(4);0])
bar(X,[0;0;0;0;Glysave(5)])
ylabel({'Mean hepatic glycogen' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
% ylim([380,460])
hold off
%%

% EGP
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'ytick',[0,0.5,1,1.5,2],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
bar(X,[EGPsave(1);0;0;0;0])
bar(X,[0;EGPsave(2);0;0;0])
bar(X,[0;0;EGPsave(3);0;0])
bar(X,[0;0;0;EGPsave(4);0])
bar(X,[0;0;0;0;EGPsave(5)])
ylabel({'Mean EGP' ; '(mg/kg/min)'},'FontSmoothing','on','fontname','Arial');
% ylim([0.85,0.96])
hold off
%%
% Gluconeogenesis
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'ytick',[0,0.5,1],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
bar(X,[Gluconeogenesissave(1);0;0;0;0])
bar(X,[0;Gluconeogenesissave(2);0;0;0])
bar(X,[0;0;Gluconeogenesissave(3);0;0])
bar(X,[0;0;0;Gluconeogenesissave(4);0])
bar(X,[0;0;0;0;Gluconeogenesissave(5)])
ylabel({'Mean gluconeogenesis' ; '(mg/kg/min)'},'FontSmoothing','on','fontname','Arial');
% ylim([0.15,0.21])
hold off

%%
%% Plot meal speeds
% Glucose AUC
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,3,6],'ytick',[2,6,10],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
plot((time-13560)/60,simFast.variablevalues(:,ismember(simFast.variables,'G'))/18,'LineWidth',LineWidthValue)
plot((time-13560)/60,simSlow.variablevalues(:,ismember(simSlow.variables,'G'))/18,'LineWidth',LineWidthValue)
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
ylabel({'Plasma glucose' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
line([0.38 0.50], [2 2],'Color','k','LineWidth',12);
xlim([0, 6])
ylim([2,10])
hold off
%%
% I
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,3,6],'ytick',[0,300,600],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
plot((time-13560)/60,simFast.variablevalues(:,ismember(simFast.variables,'I')),'LineWidth',LineWidthValue)
plot((time-13560)/60,simSlow.variablevalues(:,ismember(simSlow.variables,'I')),'LineWidth',LineWidthValue)
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
ylabel({'Plasma insulin' ; '(pM)'},'FontSmoothing','on','fontname','Arial');
line([0.38 0.50], [0 0],'Color','k','LineWidth',12);
xlim([0, 6])
ylim([0,600])
hold off
%%
%Gly
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,3,6],'ytick',[0,100,200,250,300,400],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
plot((time-13560)/60,simFast.variablevalues(:,ismember(simFast.variables,'Glycogen_liver')),'LineWidth',LineWidthValue)
plot((time-13560)/60,simSlow.variablevalues(:,ismember(simSlow.variables,'Glycogen_liver')),'LineWidth',LineWidthValue)
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
ylabel({'Hepatic glycogen' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
line([0.38 0.50], [200 200],'Color','k','LineWidth',12);
xlim([0, 6])
ylim([200,300])
hold off
%%
% EGP
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,3,6],'ytick',[0,0.5,1,1.5,2],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
plot((time-13560)/60,simFast.reactionvalues(:,ismember(simFast.reactions,'EGP')),'LineWidth',LineWidthValue)
plot((time-13560)/60,simSlow.reactionvalues(:,ismember(simSlow.reactions,'EGP')),'LineWidth',LineWidthValue)
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
ylabel({'EGP' ; '(mg/kg/min)'},'FontSmoothing','on','fontname','Arial');
line([0.38 0.50], [0.5 0.5],'Color','k','LineWidth',12);
xlim([0, 6])
ylim([0.5,2])
hold off
%%
% Gluconeogenesis
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,3,6],'ytick',[0,0.5,0.8],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
plot((time-13560)/60,simFast.reactionvalues(:,ismember(simFast.reactions,'Gluconeogenesis')),'LineWidth',LineWidthValue)
plot((time-13560)/60,simSlow.reactionvalues(:,ismember(simSlow.reactions,'Gluconeogenesis')),'LineWidth',LineWidthValue)
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
ylabel({'Gluconeogenesis' ; '(mg/kg/min)'},'FontSmoothing','on','fontname','Arial');
line([0.38 0.50], [0.2 0.2],'Color','k','LineWidth',12);
xlim([0, 6])
ylim([0.2,0.8])
hold off

%% Plot Big vs Small

X = categorical({'Obese, 110kg','Normal, 80kg'});
X = reordercats(X,{'Obese, 110kg','Normal, 80kg'});

% Glucose AUC
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
set(gca,'ytick',[0,3,6],'FontSize', 65,'fontname','Arial')%,'FontSmoothing','on')
bar(X,[Gsave(8);0])
bar(X,[0;Gsave(9)])
ylabel({'Plasma glucose' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
% ylim([3,5.7])
hold off
%%
% I
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,1,2,3,5,6,7,8,9,10],'ytick',[0,100,200,400,600,800],'FontSize', 65,'fontname','Arial')%,'FontSmoothing','on')
bar(X,[Isave(8);0])
bar(X,[0;Isave(9)])
ylabel({'Plasma insulin' ; '(pM)'},'FontSmoothing','on','fontname','Arial');
% ylim([80,240])
hold off
%%
% Gly
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,1,2,3,5,6,7,8,9,10],'ytick',[0,100,200,300,400],'FontSize', 65,'fontname','Arial')%,'FontSmoothing','on')
bar(X,[Glysave(8);0])
bar(X,[0;Glysave(9)])
ylabel({'Hepatic glycogen' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
% ylim([380,460])
hold off
%%
% EGP
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'ytick',[0,0.5,1,1.5,2],'FontSize', 65,'fontname','Arial')%,'FontSmoothing','on')
bar(X,[EGPsave(8);0])
bar(X,[0;EGPsave(9)])
ylabel({'EGP' ; '(mg/kg/min)'},'FontSmoothing','on','fontname','Arial');
% ylim([0.85,0.96])
hold off
%%
% Gluconeogenesis
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'ytick',[0,0.5,1,1.5,2],'FontSize', 65,'fontname','Arial')%,'FontSmoothing','on')
bar(X,[Gluconeogenesissave(8);0])
bar(X,[0;Gluconeogenesissave(9)])
ylabel({'Gluconeogenesis' ; '(mg/kg/min)'},'FontSmoothing','on','fontname','Arial');
ylim([0,1])
hold off