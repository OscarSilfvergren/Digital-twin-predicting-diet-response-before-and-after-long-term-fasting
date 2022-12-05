%% Create Mex och Assign params
modelName = 'ModelSimpleMetabolismFast';
model = IQMmodel([modelName '.txt']);
IQMmakeMEXmodel(model);
model = str2func(modelName);

time = linspace(0,11500,11501);
time2 = [time, fliplr(time)];
disp (' ')
disp('--- The fit is visually studied ---')
disp(' ')

%% P1 Fasting intervention
load('Silfvergren2021_ParamHealthyCalibratedP1');

[row column] = size(Silfvergren2021_ParamHealthyCalibratedp1);
clear optimizedParamTemp2

body_information  = [0 , 1, 179, 87 ];  % female, male, height, weight
meal_information  = [1, 0, 0, 0];
FedState = 1;

for i = 1:row
    optimizedParamTemp  = Silfvergren2021_ParamHealthyCalibratedp1(i,1:(column-1));
    cost                = Silfvergren2021_costfunction(Silfvergren2021_data.p1_glucoseCalibrated,Silfvergren2021_data.time,FedState,time,log(optimizedParamTemp),modelName,body_information,meal_information,190);
    optimizedParamTemp(column) = cost;
    optimizedParamTemp2(i,:)   = optimizedParamTemp;
end
Silfvergren2021_ParamHealthyCalibratedp1 = sortrows(optimizedParamTemp2,column);

for i = 1:row
    optimizedParamTemp  = AssignParameter(Silfvergren2021_ParamHealthyCalibratedp1(i,1:(column-1)), body_information, meal_information);
    
    try
        sim = model(time,[],optimizedParamTemp);
    catch error
    end
    
    if i == 1
        simBest = sim;
    end
    
    if i == 1
        maxG2 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        minG2 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
    else
        maxG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        minG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        maxG2 = max(maxG2,maxG1);
        minG2 = min(minG2,minG1);
    end
    
end

MaxP1 = maxG2;
MinP1 = minG2;

figure('Name', "Glucose fasting participant 1", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,1,2],'ytick',[0,2,3,4,5,6,7],'FontSize', 70,'fontname','Arial')
yP1    = [MaxP1', fliplr(MinP1')];
fill(time2/1440-5.33,yP1,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/1440-5.33, simBest.variablevalues(:,ismember(simBest.variables,'G'))/18,'b-.','LineWidth',LineWidthValue);
plot(Silfvergren2021_data.time(1:3)/1440-5.33,Silfvergren2021_data.p1_glucoseCalibrated(1:3),'kx','MarkerSize',80,'LineWidth',10)
plot(Silfvergren2021_data.time(1:190)/1440-5.33,Silfvergren2021_data.p1_glucoseCalibrated(1:190),'k.','MarkerSize',80)
ylabel({'Plasma glucose' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (days)",'FontSmoothing','on','fontname','Arial');
xlim([0 2])
ylim([3 6])
hold off

%% P2
load('Silfvergren2021_ParamHealthyCalibratedP2');

[row column] = size(Silfvergren2021_ParamHealthyCalibratedp2);
clear optimizedParamTemp2

% P2
body_information = [0 , 1, 178, 80 ];  % female, male, height, weight
meal_information = [1, 0, 0, 0];
FedState = 1;

for i = 1:row
    optimizedParamTemp  = Silfvergren2021_ParamHealthyCalibratedp2(i,1:(column-1));
    cost                = Silfvergren2021_costfunction(Silfvergren2021_data.p2_glucose,Silfvergren2021_data.time,FedState,time,log(optimizedParamTemp),modelName,body_information,meal_information,190);
    optimizedParamTemp(column) = cost;
    optimizedParamTemp2(i,:)   = optimizedParamTemp;
end
Silfvergren2021_ParamHealthyCalibratedp2 = sortrows(optimizedParamTemp2,column);

for i = 1:row
        optimizedParamTemp  = AssignParameter(Silfvergren2021_ParamHealthyCalibratedp2(i,1:(column-1)), body_information, meal_information);
    
    try
        sim = model(time,[],optimizedParamTemp);
    catch error
    end
    
    if i == 1
        simBest = sim;
    end
    
    if i == 1
        maxG2 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        minG2 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
    else
        maxG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        minG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        maxG2 = max(maxG2,maxG1);
        minG2 = min(minG2,minG1);
    end
    
end
MaxP2 = maxG2;
MinP2 = minG2;

%
figure('Name', "Glucose fasting participant 2", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,1,2],'ytick',[0,2,3,4,5,6,7],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
yP2    = [MaxP2', fliplr(MinP2')];
fill(time2/1440-5.33,yP2,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/1440-5.33, simBest.variablevalues(:,ismember(simBest.variables,'G'))/18,'b-.','LineWidth',LineWidthValue);
plot(Silfvergren2021_data.time(2:4)/1440-5.33,Silfvergren2021_data.p2_glucose(2:4),'kx','MarkerSize',80,'LineWidth',10)
plot(Silfvergren2021_data.time(2:190)/1440-5.33,Silfvergren2021_data.p2_glucose(2:190),'k.','MarkerSize',80)
ylabel({'Plasma glucose' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (days)",'FontSmoothing','on','fontname','Arial');
xlim([0 2])
ylim([3 6])
hold off

%% Meals
%% fed p1
load('Silfvergren2021_ParamP1fedCalibrated');

[row column] = size(Silfvergren2021_ParamP1fedCalibrated);
clear optimizedParamTemp2

% P1
body_information  = [0 , 1, 179, 87 ];  % female, male, height, weight
meal_information = [3, 0, 2.6, 25.55];
FedState = 1;

for i = 1:row
    optimizedParamTemp  = Silfvergren2021_ParamP1fedCalibrated(i,1:(column-1));
    cost                = Silfvergren2021_costfunction(Silfvergren2021_data.value_fed_p1(1:34),Silfvergren2021_data.time_fed_p1(1:34),FedState,time,log(optimizedParamTemp),modelName,body_information,meal_information,34);
    optimizedParamTemp(column) = cost;
    optimizedParamTemp2(i,:)   = optimizedParamTemp;
end
Silfvergren2021_ParamP1fedCalibrated = sortrows(optimizedParamTemp2,column);

for i = 1:row
    optimizedParamTemp  = AssignParameter(Silfvergren2021_ParamP1fedCalibrated(i,1:(column-1)), body_information, meal_information);
    
    sim = model(time,[],optimizedParamTemp);
    
    if i == 1
        simBest = sim;
    end
    
    if i == 1
        maxG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        minG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        maxG2 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        minG2 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
    else
        maxG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        minG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        maxG2 = max(maxG2,maxG1);
        minG2 = min(minG2,minG1);
    end
    
end
MaxP1 = maxG2;
MinP1 = minG2;


% Plot p1
figure('Name', "Glucose fed participant 1", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,2,4,5,6],'ytick',[0,1,2,3,4,5,6,7,8,10],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
yP1 = [MaxP1', fliplr(MinP1')];
fill((time2-7680)/60,yP1,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot((time-7680)/60, simBest.variablevalues(:,ismember(simBest.variables,'G'))/18,'b-.','LineWidth',LineWidthValue);
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
ylabel({'Plasma glucose' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
plot((Silfvergren2021_data.time_fed_p1(1:2)-7680)/60,Silfvergren2021_data.value_fed_p1(1:2),'kx','MarkerSize',80,'LineWidth',10)
plot((Silfvergren2021_data.time_fed_p1-7680)/60,Silfvergren2021_data.value_fed_p1,'k.','MarkerSize',80)
line([0 0.05], [3 3],'Color','k','LineWidth',12);
xlim([-0.1 4])
ylim([3 6])
hold off

%% fed p2
load('Silfvergren2021_ParamP2fedCalibrated');

[row column] = size(Silfvergren2021_ParamP2fedCalibrated);
clear optimizedParamTemp2

% P2
body_information = [0 , 1, 178, 80 ];  % female, male, height, weight
meal_information = [3, 0, 2.6, 25.55];
FedState = 1;

for i = 1:row
    optimizedParamTemp  = Silfvergren2021_ParamP2fedCalibrated(i,1:(column-1));
    cost                = Silfvergren2021_costfunction(Silfvergren2021_data.value_fed_p2(1:21),Silfvergren2021_data.time_fed_p2(1:21),FedState,time,log(optimizedParamTemp),modelName,body_information,meal_information,21);
    optimizedParamTemp(column) = cost;
    optimizedParamTemp2(i,:)   = optimizedParamTemp;
end
Silfvergren2021_ParamP2fedCalibrated = sortrows(optimizedParamTemp2,column);

for i = 1:row
        optimizedParamTemp  = AssignParameter(Silfvergren2021_ParamP2fedCalibrated(i,1:(column-1)), body_information, meal_information);
    
    try
        sim = model(time,[],optimizedParamTemp);
    catch error
    end
    
    if i == 1
        simBest = sim;
    end
    
    if i == 1
        maxG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        minG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        maxG2 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        minG2 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
    else
        maxG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        minG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        maxG2 = max(maxG2,maxG1);
        minG2 = min(minG2,minG1);
    end
    
end
MaxP1 = maxG2;
MinP1 = minG2;

% Plot p2
figure('Name', "Glucose fed participant 2", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,2,4,5,6],'ytick',[0,1,2,3,4,5,6,7,8,10],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
yP1    = [MaxP1', fliplr(MinP1')];
fill((time2-7680)/60,yP1,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot((time-7680)/60, simBest.variablevalues(:,ismember(simBest.variables,'G'))/18,'b-.','LineWidth',LineWidthValue);
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
ylabel({'Plasma glucose' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
plot((Silfvergren2021_data.time_fed_p2(1:2)-7680)/60,Silfvergren2021_data.value_fed_p2(1:2),'kx','MarkerSize',80,'LineWidth',10)
plot((Silfvergren2021_data.time_fed_p2-7680)/60,Silfvergren2021_data.value_fed_p2,'k.','MarkerSize',80)
line([0 0.05], [3 3],'Color','k','LineWidth',12);
xlim([-0.1 4])
ylim([3, 6])
hold off

%% unfed p1
load('Silfvergren2021_ParamP1unfedCalibrated');

[row column] = size(Silfvergren2021_ParamP1unfedCalibrated);
clear optimizedParamTemp2

% P1
body_information  = [0 , 1, 179, 87 ];  % female, male, height, weight
meal_information  = [3, 0, 2.6, 25.55];

FedState = 0;

for i = 1:row
    optimizedParamTemp  = Silfvergren2021_ParamP1unfedCalibrated(i,1:(column-1));
    cost                = Silfvergren2021_costfunction(Silfvergren2021_data.Value_fasted_p1(1:56),Silfvergren2021_data.time_fastedStart_p1(1:56),FedState,time,log(optimizedParamTemp),modelName,body_information,meal_information,56);
    optimizedParamTemp(column) = cost;
    optimizedParamTemp2(i,:)   = optimizedParamTemp;
end
Silfvergren2021_ParamP1unfedCalibrated = sortrows(optimizedParamTemp2,column);

for i = 1:row
    optimizedParamTemp  = AssignParameter(Silfvergren2021_ParamP1unfedCalibrated(i,1:(column-1)), body_information, meal_information);
    
    try
        sim = model(time,[],optimizedParamTemp);
    catch error
    end
    
    if i == 1
        simBest = sim;
    end
    
    if i == 1
        maxG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        minG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        maxG2 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        minG2 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
    else
        maxG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        minG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        maxG2 = max(maxG2,maxG1);
        minG2 = min(minG2,minG1);
    end
    
end
MaxP1 = maxG2;
MinP1 = minG2;

% Plot p1
figure('Name', "Glucose unfed participant 1", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,2,4,5,6],'ytick',[2,4,6],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
yP1    = [MaxP1', fliplr(MinP1')];
fill((time2-10560)/60,yP1,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot((time-10560)/60, simBest.variablevalues(:,ismember(simBest.variables,'G'))/18,'b-.','LineWidth',LineWidthValue);
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
ylabel({'Plasma glucose' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
plot((Silfvergren2021_data.time_fastedStart_p1(1:2)-10560)/60,Silfvergren2021_data.Value_fasted_p1(1:2),'kx','MarkerSize',80,'LineWidth',10)
plot((Silfvergren2021_data.time_fastedStart_p1-10560)/60,Silfvergren2021_data.Value_fasted_p1,'k.','MarkerSize',80)
line([0.2 0.25], [2 2],'Color','k','LineWidth',12);
xlim([-0.1 4])
ylim([2, 6])
hold off


%% unfed p2
load('Silfvergren2021_ParamP2unfedCalibrated');

[row column] = size(Silfvergren2021_ParamP2unfedCalibrated);
clear optimizedParamTemp2

% P2
body_information = [0 , 1, 178, 80 ];  % female, male, height, weight
meal_information = [3, 0, 2.6, 25.55];
FedState = 0;

for i = 1:row
    optimizedParamTemp  = Silfvergren2021_ParamP2unfedCalibrated(i,1:(column-1));
    cost                = Silfvergren2021_costfunction(Silfvergren2021_data.Value_fasted_p2(1:44),Silfvergren2021_data.time_fastedStart_p2(1:44),FedState,time,log(optimizedParamTemp),modelName,body_information,meal_information,44);
    optimizedParamTemp(column) = cost;
    optimizedParamTemp2(i,:)   = optimizedParamTemp;
end
Silfvergren2021_ParamP2unfedCalibrated = sortrows(optimizedParamTemp2,column);

for i = 1:row
    optimizedParamTemp  = AssignParameter(Silfvergren2021_ParamP2unfedCalibrated(i,1:(column-1)), body_information, meal_information);
    
    try
    sim = model(time,[],optimizedParamTemp);
    end
    
    if i == 1
        simBest = sim;
    end
    
    if i == 1
        maxG2 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        minG2 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
    else
        maxG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        minG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        maxG2 = max(maxG2,maxG1);
        minG2 = min(minG2,minG1);
    end
    
    if i == 1
        maxGNG2 = sim.reactionvalues(:,ismember(sim.reactions,'Gluconeogenesis'));
        minGNG2 = sim.reactionvalues(:,ismember(sim.reactions,'Gluconeogenesis'));
    else
        maxGNG1 = sim.reactionvalues(:,ismember(sim.reactions,'Gluconeogenesis'));
        minGNG1 = sim.reactionvalues(:,ismember(sim.reactions,'Gluconeogenesis'));
        maxGNG2 = max(maxGNG2,maxGNG1);
        minGNG2 = min(minGNG2,minGNG1);
    end
    
    if i == 1
        maxPT2 = sim.reactionvalues(:,ismember(sim.reactions,'PyruvateTranslocase'));
        minPT2 = sim.reactionvalues(:,ismember(sim.reactions,'PyruvateTranslocase'));
    else
        maxPT1 = sim.reactionvalues(:,ismember(sim.reactions,'PyruvateTranslocase'));
        minPT1 = sim.reactionvalues(:,ismember(sim.reactions,'PyruvateTranslocase'));
        maxPT2 = max(maxPT2,maxPT1);
        minPT2 = min(minPT2,minPT1);
    end
    
    if i == 1
        maxRatio2 = sim.reactionvalues(:,ismember(sim.reactions,'aaIntoLiver_Meal'))./sim.reactionvalues(:,ismember(sim.reactions,'aaTransportation'));
        minRatio2 = sim.reactionvalues(:,ismember(sim.reactions,'aaIntoLiver_Meal'))./sim.reactionvalues(:,ismember(sim.reactions,'aaTransportation'));
    else
        maxRatio1 = sim.reactionvalues(:,ismember(sim.reactions,'aaIntoLiver_Meal'))./sim.reactionvalues(:,ismember(sim.reactions,'aaTransportation'));
        minRatio1 = sim.reactionvalues(:,ismember(sim.reactions,'aaIntoLiver_Meal'))./sim.reactionvalues(:,ismember(sim.reactions,'aaTransportation'));
        maxRatio2 = max(maxRatio2,maxRatio1);
        minRatio2 = min(minRatio2,minRatio1);
    end
    
    if i == 1
        maxProteinDigestion2 = sim.reactionvalues(:,ismember(sim.reactions,'ProteinDigestion'));
        minProteinDigestion2 = sim.reactionvalues(:,ismember(sim.reactions,'ProteinDigestion'));
    else
        maxProteinDigestion1 = sim.reactionvalues(:,ismember(sim.reactions,'ProteinDigestion'));
        minProteinDigestion1 = sim.reactionvalues(:,ismember(sim.reactions,'ProteinDigestion'));
        maxProteinDigestion2 = max(maxProteinDigestion2,maxProteinDigestion1);
        minProteinDigestion2 = min(minProteinDigestion2,minProteinDigestion1);
    end
    
    
end

MaxP1 = maxG2;
MinP1 = minG2;

% Plot p2
figure('Name', "Glucose unfed participant 2", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,2,4,5,6],'ytick',[2,4,6],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
yP1    = [MaxP1', fliplr(MinP1')];
fill((time2-10560)/60,yP1,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot((time-10560)/60, simBest.variablevalues(:,ismember(simBest.variables,'G'))/18,'b-.','LineWidth',LineWidthValue);
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
ylabel({'Plasma glucose' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
plot((Silfvergren2021_data.time_fastedStart_p2(1:2)-10560)/60,Silfvergren2021_data.Value_fasted_p2(1:2),'kx','MarkerSize',80,'LineWidth',10)
plot((Silfvergren2021_data.time_fastedStart_p2-10560)/60,Silfvergren2021_data.Value_fasted_p2,'k.','MarkerSize',80)
line([0.2 0.25], [2 2],'Color','k','LineWidth',12);
xlim([-0.1 4])
ylim([2, 6])
hold off

figure('Name', "GNG unfed p2", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,1,2],'ytick',[0.2,0.6,1],'FontSize', 65,'fontname','Arial')
yP1    = [minGNG2', fliplr(maxGNG2')];
fill(time2/1440-5.33,yP1,'b','FaceAlpha',0.2,'EdgeAlpha',0);
line([0 0.05], [0 0],'Color','k','LineWidth',12);
line([2 2.05], [0 0],'Color','k','LineWidth',12);
plot(time/1440-5.33, simBest.reactionvalues(:,ismember(simBest.reactions,'Gluconeogenesis')),'b-.','LineWidth',LineWidthValue);
ylabel({'Gluconeogenesis' ; '(mg/kg/min)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (days)",'FontSmoothing','on','fontname','Arial');
xlim([-0.1 2.5])
ylim([0 1])
hold off

figure('Name', "PyruvateTranslocase", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,1,2],'ytick',[8,9,10,12],'FontSize', 65,'fontname','Arial')
yP1    = [maxPT2', fliplr(minPT2')];
fill(time2/1440-5.33,yP1*1000,'b','FaceAlpha',0.2,'EdgeAlpha',0);
line([0 0.05], [8 8],'Color','k','LineWidth',12);
line([2 2.05], [8 8],'Color','k','LineWidth',12);
plot(time/1440-5.33, simBest.reactionvalues(:,ismember(simBest.reactions,'PyruvateTranslocase'))*1000,'b-.','LineWidth',LineWidthValue);
ylabel({'Pyruvate from body' ; 'into liver (ug/kg/min)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (days)",'FontSmoothing','on','fontname','Arial');
xlim([-0.1 2.5])
ylim([8 12])
hold off

figure('Name', "digested AA ratio", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,1,2],'ytick',[0,50,100],'FontSize', 55,'fontname','Arial')
yP1    = [minRatio2', fliplr(maxRatio2')];
fill(time2(2:end-1)/1440-5.33,yP1(2:end-1)*100,'b','FaceAlpha',0.2,'EdgeAlpha',0);
bestRatio = sim.reactionvalues(:,ismember(sim.reactions,'aaIntoLiver_Meal'))./sim.reactionvalues(:,ismember(sim.reactions,'aaTransportation'));
plot(time(2:end-1)/1440-5.33, bestRatio(2:end-1)*100,'b-.','LineWidth',LineWidthValue);
line([0 0.05], [0 0],'Color','k','LineWidth',12);
line([2 2.05], [0 0],'Color','k','LineWidth',12);
ylabel({'Amino acids catabolized' ; 'into TCA or pyruvate (%)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (days)",'FontSmoothing','on','fontname','Arial');
xlim([-0.1 2.5])
ylim([0 100])
hold off

figure('Name', "Digestion of protein into AA", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,1,2],'ytick',[0,1.5,3],'FontSize', 50,'fontname','Arial')
yP1    = [maxProteinDigestion2', fliplr(minProteinDigestion2')];
fill(time2/1440-5.33,yP1,'b','FaceAlpha',0.2,'EdgeAlpha',0);
line([0 0.05], [0 0],'Color','k','LineWidth',12);
line([2 2.05], [0 0],'Color','k','LineWidth',12);
plot(time/1440-5.33, sim.reactionvalues(:,ismember(sim.reactions,'ProteinDigestion')),'b-.','LineWidth',LineWidthValue);
ylabel({'Release of amino acids from' ; 'ingested protein(mg/kg/min)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (days)",'FontSmoothing','on','fontname','Arial');
xlim([-0.1 2.5])
ylim([0 3])
hold off