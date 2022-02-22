%% Simulate

% old model
modelName = 'DallamanLerche';
model = IQMmodel([modelName '.txt']);
IQMmakeMEXmodel(model);
model = str2func(modelName);
load('Dallaman_paramsHealthy');
simDallamanHealthy = model(time,[],Dallaman_paramsHealthy);

% new model
modelName = 'ModelSimpleMetabolismSteadyState';
model = IQMmodel([modelName '.txt']);
IQMmakeMEXmodel(model);
model = str2func(modelName);
load('Lerche2009_Param');
simLerche2009Healthy = model(time,[], Lerche2009_Param);

%% Plot fast

%Glucose
figure('Name', "Compare Glucose fasting", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
set(gca,'xtick',[0,1,2],'ytick',[3,4,5],'FontSize', 70,'FontSmoothing','on','fontname','Arial')
xlabel("Time (days)",'FontSmoothing','on');
ylabel({'Plasma glucose' ; '(mM)'},'FontSmoothing','on');
%ylim([0 250])
plot(time/1440-4.98,simLerche2009Healthy.variablevalues(:,ismember(simLerche2009Healthy.variables,'G'))/18,'Color',colorGrey,'LineWidth',LineWidthValue)
plot(time/1440-4.98,simDallamanHealthy.reactionvalues(:,ismember(simDallamanHealthy.reactions,'G'))/18,'LineWidth',LineWidthValue,'Color',colorOrange)
errorbar(Lerche2009_data.time_glucose_healthy/1440-4.98,Lerche2009_data.glucose_healthy/18,Lerche2009_data.glucoseSEM_healthy/18, " k . ",'MarkerSize',8,'LineWidth',4);
ylim([3 5.2])
xlim([0 2])
hold off

% Insulin
figure('Name', "Compare Insulin fasting", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
set(gca,'xtick',[0,1,2],'ytick',[0,20,40],'FontSize', 70,'FontSmoothing','on','fontname','Arial')
xlabel("Time (days)",'FontSmoothing','on');
ylabel({'Plasma insulin' ; '(pM)'},'FontSmoothing','on');
xlim([100 180])
plot(time/1440-4.98,simLerche2009Healthy.variablevalues(:,ismember(simLerche2009Healthy.variables,'I')),'Color',colorGrey,'LineWidth',LineWidthValue)
plot(time/1440-4.98,simDallamanHealthy.reactionvalues(:,ismember(simDallamanHealthy.reactions,'I')),'LineWidth',LineWidthValue,'Color',colorOrange)
errorbar(Lerche2009_data.time_insulin_healthy/1440-4.98,Lerche2009_data.insulin_healthy,Lerche2009_data.insulinSEMFIXED_healthy, " k . ",'MarkerSize',8,'LineWidth',4);
ylim([0 40])
xlim([0 2])
hold off